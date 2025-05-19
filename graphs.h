#ifndef GRAPH_H_
#define GRAPH_H_

#include <stdio.h>
#include <vector>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <unordered_map>
#include <map>
#include <cctype>
#include <algorithm>
#include <random>
#include <set>
#include <queue>
/*#include <thread>
#include <atomic>
#include <mutex>*/

inline std::vector<double> sample_vector_C(const std::vector<double>& input, int N, bool randomSample) {
    std::vector<double> result;

    if (input.size() <= N) {
        return input; // Return whole input if it's already small enough
    }

    if (randomSample) {
        std::vector<size_t> indices(input.size());
        std::iota(indices.begin(), indices.end(), 0); // Fill with 0..n-1
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(indices.begin(), indices.end(), gen);

        for (size_t i = 0; i < N; ++i) {
            result.push_back(input[indices[i]]);
        }
    } else {
        result.insert(result.end(), input.begin(), input.begin() + N);
    }

    return result;
}

class InstanceData {
public:

    std::unordered_map<std::string, std::vector<std::string>> string_columns;
    std::unordered_map<std::string, std::vector<double>> double_columns;

    std::string get_string(const std::string& key, size_t index) const {
        return string_columns.at(key).at(index);
    }

    std::vector<std::string> get_string_array(const std::string& key) const {
        return string_columns.at(key);
    }

    double get_double(const std::string& key, size_t index) const {
        return double_columns.at(key).at(index);  
    }

    std::vector<double> get_double_array(const std::string& key) const {
        if (double_columns.find(key) != double_columns.end()) {
            return double_columns.at(key);
        }
        else {
            std::cerr << "Column " << key << " does not exist." << std::endl;
            return {};
        }
    }

    void add_column(const std::string& name, const std::vector<std::string>& column) {
        string_columns[name] = column;
    }

    void add_column(const std::string& name, const std::vector<double>& column) {
        double_columns[name] = column;
    }

};

struct edges_matrix {
    double capacity;
    double flow;

    // Constructor
    edges_matrix() : capacity(0), flow(0) {}
    edges_matrix(double c, double fl) : capacity(c), flow(fl) {}
};

struct NameMapHelp {
    int to;
    int ptr;

    NameMapHelp() : to(0), ptr(0) {}
    NameMapHelp(int t, int p) : to(t), ptr(p) {}
};

class Graph {
    public:
    void init_focs(InstanceData instance, int timeStep, int instancesize, bool randomize) {

        std::vector<double> average_power_W;
        std::vector<double> t0_x;
        std::vector<double> t1_x;
        std::vector<double> total_energy_Wh;
        std::vector<double> maxPower;

        std::vector<double> jobs_cap;
        std::vector<double> jobs_demand;
        std::vector<int> jobs_departure;
        std::vector<int> jobs_arrival;
        std::vector<int> breakpoints;

        int n;

        timestep = timeStep;
        //init all

        average_power_W = sample_vector_C(instance.get_double_array("average_power_W"), instancesize, randomize);
        total_energy_Wh = sample_vector_C(instance.get_double_array("total_energy_Wh"), instancesize, randomize);
        maxPower = sample_vector_C(instance.get_double_array("maxPower"), instancesize, randomize);

        char t0_name[10], t1_name[10];
        snprintf(t0_name, sizeof(t0_name), "t0_%d", timestep);
        snprintf(t1_name, sizeof(t1_name), "t1_%d", timestep);

        t0_x = sample_vector_C(instance.get_double_array(t0_name), instancesize, randomize);
        t1_x = sample_vector_C(instance.get_double_array(t1_name), instancesize, randomize);
        
        n = average_power_W.size();
        for (int j = 0; j < n; ++j) {
            jobs.push_back(j);
            if (average_power_W[j] > maxPower[j])
                jobs_cap.push_back(22);
            else
                jobs_cap.push_back(maxPower[j] / 1000.0);
            jobs_demand.push_back(total_energy_Wh[j] / 1000.0);
            jobs_departure.push_back(t0_x[j]);
            jobs_arrival.push_back(t1_x[j]);
        }

        std::set<int> unique_breakpoints(t0_x.begin(), t0_x.end());
        unique_breakpoints.insert(t1_x.begin(), t1_x.end());

        breakpoints.assign(unique_breakpoints.begin(), unique_breakpoints.end());

        for (size_t i = 0; i < breakpoints.size() - 1; i++) {
            intervals_start.push_back(breakpoints[i]);
            intervals_end.push_back(breakpoints[i+1]);
            I_a.push_back(i);
            len_i.push_back((breakpoints[i+1] - breakpoints[i]) * timestep);
        }

        find_J();
        find_J_inverse(t0_x, t1_x);

        //HERE INSTANTIATION OF MATRIX

        iX = jobs.size() + 1;
        tX = iX + I_a.size();
        NameMap.resize(tX);
        FastNameMap.resize(tX);

        //First add all names to NameMap and ReverseNameMap
        ReverseNameMapF.push_back(0);
        ReverseNameMapS.push_back(0);
        int edgeIndex = 1;      //Starting at 1 so if it defaults to 0 we know it couldnt find the edge.

        // Add all edge names to NameMap and ReverseNameMap
        for (int j : jobs) {
            int from = sX;
            int to = jX + j;
            FastNameMap[from].push_back(NameMapHelp(to, edgeIndex));
            NameMap[from][to] = edgeIndex;
            ReverseNameMapF.push_back(from);
            ReverseNameMapS.push_back(to);
            edgeIndex++;
        }

        for (int j : jobs) {
            int from = jX + j;
            double cap_j = jobs_cap[j];
            
            for (int i : J_inverse[from]) {
                int to = iX + i;
                FastNameMap[from].push_back(NameMapHelp(to, edgeIndex));
                NameMap[from][to] = edgeIndex;
                ReverseNameMapF.push_back(from);
                ReverseNameMapS.push_back(to);
                edgeIndex++;
            }
        }

        for (int i : I_a) {
            int from = iX + i;
            int to = tX;
            FastNameMap[from].push_back(NameMapHelp(to, edgeIndex));
            NameMap[from][to] = edgeIndex;
            ReverseNameMapF.push_back(from);
            ReverseNameMapS.push_back(to);
            edgeIndex++;
        }


        //dedicate enough space to all matrices.
        G = std::vector<edges_matrix>(edgeIndex);
        G_rk = G;
        f_matrix = G;
        //now we can access each matrix via: G[NameMap[sX]["j0"]].capacity = 1;


        //focs_instance_to_network
        for (int j : jobs) {
            int from = sX;
            int to = jX + j;
            double capacity = jobs_demand[j];
            G[NameMap[from][to]].capacity = capacity;
        }
        
        for (int j : jobs) {
            int from = jX + j;
            double cap_j = jobs_cap[j];
        
            for (int i : J_inverse[from]) {
                int to = iX + i;
                double capacity = cap_j * len_i[i] / timeBase;
                G[NameMap[from][to]].capacity = capacity;
            }
        }
        
        for (int i : I_a) {
            int from = iX + i;
            int to = tX;
            double capacity = 0.0;  // If no value is given, use 0 or determine based on context
            G[NameMap[from][to]].capacity = capacity;
        }     

        //std::cout << "TimeComplexity is: O(" << std::max(pow(jobs.size(),4), pow(I_a.size(),4))  << ") \n";
        //std::cout << "TimeComplexity with Bi-directional BFS: O(" << std::max(pow(jobs.size(),2), pow(I_a.size(),2))  << ") \n";


        build_reverse_adj();

        G_r = G;
        total_demand_r = Get_M(G_r);  
        selfterminate = false;
        it = 0;
    }

    double totalEDMONDKARP = 0;
    void solve_focs() {
        reset_caps(f_matrix);
        reset_flows(f_matrix);
        while (!selfterminate) {
            if(it == 0) {
                G_rk = G_r;
            }
            update_network_capacities_g();
            //Edmonds_Karp();
            clock_t z1 = clock();
            Max_flow_solver();
            clock_t z2 = clock();
            totalEDMONDKARP += (double)(z2-z1);
            MaxDiff = Get_M(G_rk) - flow_val; 
            
            if (total_demand_r-flow_val < err) {
                //end round

                //Check for subcrit. For some instances, there are still active intervals that only now don't reach the max anymore
                subCrit_mask.clear();
                int toKey = tX;
                for (int i : I_a) {
                    int fromKey = iX + i;
                    edges_matrix TheEdge = G_rk[NameMap[fromKey][toKey]];

                    bool isCritical = (TheEdge.capacity - TheEdge.flow > err);
                    subCrit_mask.push_back(isCritical);
                }


                subCrit.clear();
                for (size_t i = 0; i < I_a.size(); ++i) {
                    if (subCrit_mask[i]) {
                        subCrit.push_back(I_a[i]);
                    }
                }

                //Update I_p
                I_p.insert(I_p.end(), subCrit.begin(), subCrit.end());

                //Update I_crit
                std::vector<int> temp;

                for (size_t i = 0; i < I_a.size(); ++i) {
                    if (!subCrit_mask[i]) {
                        temp.push_back(I_a[i]);
                    }
                }
                std::sort(temp.begin(), temp.end());
                I_crit.push_back(temp);
                if (I_p.size() == 0) {
                    selfterminate = true;
                    Add_to_f_end(G_rk);
                }
                else {
                    I_crit_r = I_crit.back();
                    AddingF = true;
                    reduce_network(I_crit_r, G_rk);
                    AddingF = false;
                    I_a = I_p;                     // Copy the contents of I_p to I_a
                    std::sort(I_a.begin(), I_a.end());  // Sort I_a in ascending order
                    I_p.clear();

                    total_demand_r = Get_M(G_r)-flow_val_saved;

                    flow_val_saved += flow_val;
                
                    rd++;
                    it = 0;
                }
            }
            else {
                
                subCrit_mask.clear();
                int toKey = tX;
                for (int i : I_a) {
                    int fromKey = iX + i;
                    edges_matrix TheEdge = G_rk[NameMap[fromKey][toKey]];
                    bool isCritical = (TheEdge.capacity - TheEdge.flow > err);
                    subCrit_mask.push_back(isCritical);
                }

                subCrit.clear();
                for (size_t i = 0; i < I_a.size(); ++i) {
                    if (subCrit_mask[i]) {
                        subCrit.push_back(I_a[i]);
                    }
                }

                reduce_network(subCrit, G_rk);

                total_demand_r = Get_M(G_rk)-flow_val_saved;

                I_p.insert(I_p.end(), subCrit.begin(), subCrit.end());

                std::vector<int> new_I_a;
                for (size_t i = 0; i < I_a.size(); ++i) {
                    if (!subCrit_mask[i]) {
                        new_I_a.push_back(I_a[i]);
                    }
                }
                std::sort(new_I_a.begin(), new_I_a.end());
                I_a = new_I_a;
                it++;
            } 
            if (it > 10) {
                selfterminate = true;
            }
        }
    }

    void objective() {
        int m = intervals_start.size();
        std::vector<double> p_i(m);
        for (int i = 0; i < m; i++) {
            int Ikey = iX + i;
            double flow = f_matrix[NameMap[Ikey][tX]].flow;
            p_i[i] = (flow / len_i[i]) * timeBase;
        }
        std::vector<double> powerSquare(m);
        for (int i = 0; i < m; ++i) {
            powerSquare[i] = (p_i[i] * p_i[i]) * (len_i[i] / timestep);
        }

        objNormalized = std::accumulate(powerSquare.begin(), powerSquare.end(), 0.0);
        std::cout << "Objective value C++ = " << objNormalized << "\n";
        //return objNormalized;
    }


    private: 

    void find_J() {
        for (int i : I_a) {
            int i_key = iX + i;
            std::vector<int> related_jobs;
        
            for (int j : jobs) {
                int j_key = jX + j;
                const auto& job_list = J_inverse[j_key];
                if (std::find(job_list.begin(), job_list.end(), i) != job_list.end()) {
                    related_jobs.push_back(j);
                }
            }
        
            J[i_key] = related_jobs;
        }        
    }

    void find_J_inverse(std::vector<double> t0_x, std::vector<double> t1_x) {
        for (int j = 0; j < jobs.size(); j++) {
            int key = jX + j;
            int t0 = t0_x[j];
            int t1 = t1_x[j];
    
            // Find index of t0 in intervals_start
            auto it_start = std::find(intervals_start.begin(), intervals_start.end(), t0);
            auto it_end = std::find(intervals_end.begin(), intervals_end.end(), t1);
    
            if (it_start != intervals_start.end() && it_end != intervals_end.end()) {
                int idx_start = std::distance(intervals_start.begin(), it_start);
                int idx_end = std::distance(intervals_end.begin(), it_end);
    
                std::vector<int> indices;
                for (int i = idx_start; i <= idx_end; ++i) {
                    indices.push_back(i);
                }
    
                J_inverse[key] = indices;
            } else {
                std::cerr << "Error: t0 or t1 not found in interval lists for job " << j << "\n";
            }
        }
    }

    void build_reverse_adj() {
    // Step 1: Find max 'to' node index
    int max_index = 0;
    for (int from = 0; from < NameMap.size(); ++from) {
        for (const auto& [to, edgeIdx] : NameMap[from]) {
            if (to >= max_index) {
                max_index = to + 1;
            }
        }
    }

    // Step 2: Count how many edges point to each node
    std::vector<int> in_degree(max_index, 0);
    for (int from = 0; from < NameMap.size(); ++from) {
        for (const auto& [to, _] : NameMap[from]) {
            ++in_degree[to];
        }
    }

    // Step 3: Initialize and reserve space for reverse_adj
    reverse_adj.clear();
    reverse_adj.resize(max_index);
    for (int i = 0; i < max_index; ++i) {
        reverse_adj[i].reserve(in_degree[i]);
    }

    // Step 4: Fill reverse_adj
    for (int from = 0; from < NameMap.size(); ++from) {
        for (const auto& [to, edgeIdx] : NameMap[from]) {
            reverse_adj[to].push_back(edgeIdx);
        }
    }
}

    void reduce_network(std::vector<int> crit_r, std::vector<edges_matrix>& graph_rK) {
        for (int i : crit_r) {
            int Ikey = iX + i;
            double remove_from_s_to_jX = 0;

            for (int j : reverse_adj[Ikey]) {   //if goes towards iX
                remove_from_s_to_jX = graph_rK[j].flow;
                graph_rK[NameMap[sX][ReverseNameMapF[j]]].capacity -= remove_from_s_to_jX;
                graph_rK[NameMap[sX][ReverseNameMapF[j]]].flow -= remove_from_s_to_jX;;
                if(AddingF) {  // if F is initiated in another round, add the flows to F.
                    f_matrix[j].flow += graph_rK[j].flow;
                }
            }

            for (const auto& [to, j] : NameMap[Ikey]) { //for goes from iX
                if(AddingF) {  // if F is initiated in another round, add the flows to F.
                    f_matrix[j].flow += graph_rK[j].flow;
                }
                graph_rK[j].capacity = 0;
                graph_rK[j].flow = 0;
            }
        }
    }

    void Add_to_f_end(std::vector<edges_matrix>& graph_rK) {
        for (int i = 0; i < graph_rK.size(); i++) {
            f_matrix[i].flow += graph_rK[i].flow;
        }
    }

    void reset_flows(std::vector<edges_matrix>& GRAPH) {
        for (int i = 0; i < GRAPH.size(); i++) {
            GRAPH[i].flow = 0;
        }
    }

    void reset_caps(std::vector<edges_matrix>& GRAPH) {
        for (int i = 0; i < GRAPH.size(); i++) {
            GRAPH[i].capacity = 0;
        }
    }

    void update_network_capacities_g() {
        double demand = 0;
        double demand_normalized = 0;
        if (it > 0) {
            if (rd == 0) {
                demand = MaxDiff;
            }
            else {
                demand = MaxDiff-flow_val_saved;
            }
            
            demand_normalized = demand / length_sum_intervals(I_a, len_i);
            for (int i : I_a) {
                int Ikey = iX + i;
                G_rk[NameMap[Ikey][tX]].capacity += demand_normalized * len_i[i];
            }
        }
        else {
            if (rd == 0) {
                demand = Get_M(G_rk);
            }
            else {
                demand = Get_M(G_rk)-flow_val_saved;
            }
            demand_normalized = demand / length_sum_intervals(I_a, len_i);
            for (int i : I_a) {
                int Ikey = iX + i;
                G_rk[NameMap[Ikey][tX]].capacity = demand_normalized * len_i[i];
            }
        }
    }
    
    void Max_flow_solver() {
        clock_t z1 = clock();
        reset_flows(G_rk);   // Reset all flows to 0 for a DFS effect
        Edmonds_Karp_Bidirectional();

        clock_t z2 = clock();
        //EDMONDSKARPTIME += (double)(z2-z1);
    }

    //Breadth First Search
    bool bfs(std::vector<int>& parent, const int& source, const int& sink, const std::vector<edges_matrix>& graph) {
        static std::deque<int> q;
        q.push_front(source);

        std::fill(parent.begin()+1, parent.end(), -1);    //fill parent with -1 for all values.
        int u = 0;
        while (!q.empty()) {
            u = q.front();
            q.pop_front();

            for (const auto& [to, edgeIdx] : FastNameMap[u]) {
                const edges_matrix& a = graph[edgeIdx];
                // `to` is the destination node name (e.g., "j0")
                // `edgeIdx` is the unique index for edge (sX -> to)
                if(a.capacity - a.flow > err){
                    if (to == sink) {
                        parent[to] = u;
                        q.clear();
                        return true;
                    }
                    else if (parent[to] == -1) {
                        q.push_back(to);   //using a stack instead of a queue, LIFO instead of FIFO, results in a +- 8% speedup
                        parent[to] = u;
                    }
                }
            }
            for (int idx : reverse_adj[u]) {
                int from = ReverseNameMapF[idx];
                if (parent[from] == -1 && graph[idx].flow > err) {
                    q.push_back(from);
                    parent[from] = ReverseNameMapS[idx];
                }
            }
        }
        return false;
    }

    bool bidirectional_bfs(std::vector<int>& parent, std::vector<int>& parent_rev, int& meet_node, const int& source, const int& sink, const std::vector<edges_matrix>& graph) {
        static std::deque<int> q_fwd, q_bwd;
        parent.assign(parent.size(), -1);
        parent_rev.assign(parent_rev.size(), -1);

        q_fwd.push_front(source);
        parent[source] = source;

        //ATOMIC BOOL NOT FOUND

        q_bwd.push_front(sink);
        parent_rev[sink] = sink;
        while (!q_fwd.empty() && !q_bwd.empty()) {

            // Expand forward search
            if (!q_fwd.empty()) {   //&&NOT FOUND
                int u = q_fwd.front();  
                q_fwd.pop_front();

                for (const auto& [to, edgeIdx] : FastNameMap[u]) {
                    const edges_matrix& a = graph[edgeIdx];
                    if (a.capacity - a.flow > err && parent[to] == -1) {    //&&NOT FOUND
                        parent[to] = u; //mtx lock
                        q_fwd.push_back(to);    //mtx lock
                                                //start new thread if pushed_back.??
                        if (parent_rev[to] != -1) {
                            meet_node = to;
                            q_fwd.clear();  //mtx lock
                            q_bwd.clear();  //mtx lock
                            return true;    
                        }
                    }
                }
                for (int idx : reverse_adj[u]) {
                    int from = ReverseNameMapF[idx];
                    if (parent[from] == -1 && graph[idx].flow > err) {  //&&NOT FOUND
                        parent[from] = ReverseNameMapS[idx];    //mtx lock
                        q_fwd.push_back(from);  //mtx lock
                        if (parent_rev[from] != -1) {
                            meet_node = from;
                            q_fwd.clear();  //mtx lock
                            q_bwd.clear();  //mtx lock
                            return true;
                        }
                    }
                }
            }

            // Expand backward search
            if (!q_bwd.empty()) {   //&&NOT FOUND
                int v = q_bwd.front();
                q_bwd.pop_front();

                for (int idx : reverse_adj[v]) {
                    int from = ReverseNameMapF[idx];
                    const edges_matrix& a = graph[idx];
                    if (a.capacity - a.flow > err && parent_rev[from] == -1) {  //&&NOT FOUND
                        parent_rev[from] = ReverseNameMapS[idx];    //mtx lock
                        q_bwd.push_back(from);  //mtx lock
                        if (parent[from] != -1) {
                            meet_node = from;
                            q_fwd.clear();  //mtx lock
                            q_bwd.clear();  //mtx lock
                            return true;
                        }
                    }
                }

                if(v != sink) {
                    for (const auto& [to, edgeIdx] : FastNameMap[v]) {
                        //if(to != sink){
                            if (graph[edgeIdx].flow > err && parent_rev[to] == -1) {  // Reverse in backward search //&&NOT FOUND
                                parent_rev[to] = v; //mtx lock
                                q_bwd.push_back(to);    //mtx lock
                                if (parent[to] != -1) {
                                    meet_node = to;
                                    q_fwd.clear();  //mtx lock
                                    q_bwd.clear();  //mtx lock
                                    return true;
                                }
                            }
                        //}
                    }
                }
            }
        }
        //THREADS JOIN()
        //IF FOUND RETURN TRUE ELSE FALSE. 
        return false;
    }

    //std::atomic<bool> FoundMeet = false;
    //std::atomic<int> MeetnodeAtomic;
    /*bool FoundMeet = false;
    int MeetnodeAtomic = -1;

    std::mutex mtx_parent;
    std::mutex mtx_parent_rev;
    std::mutex mtx_q_fwd;
    std::mutex mtx_q_bwd;

    bool bidirectional_bfs_threaded(std::vector<int>& parent, std::vector<int>& parent_rev, int& meet_node, const int& source, const int& sink, const std::vector<edges_matrix>& graph) {
        std::deque<int> q_fwd, q_bwd;
        parent.assign(parent.size(), -1);
        parent_rev.assign(parent_rev.size(), -1);

        q_fwd.push_front(source);
        parent[source] = source;

        q_bwd.push_front(sink);
        parent_rev[sink] = sink;
        //MeetnodeAtomic.store(-1);
        //FoundMeet.store(false);
        bool found = FoundMeet; //FoundMeet.load();
        while (!q_fwd.empty() && !q_bwd.empty() && !found) {
            forward_search(parent, parent_rev, q_fwd, q_bwd, graph);
            backward_search(parent, parent_rev, q_fwd, q_bwd, sink, graph);
            //std::thread t1(&Graph::forward_search, this, std::ref(parent), std::ref(parent_rev), std::ref(q_fwd), std::ref(q_bwd), std::ref(graph));
            //std::thread t2(&Graph::backward_search, this, std::ref(parent), std::ref(parent_rev), std::ref(q_fwd), std::ref(q_bwd), std::ref(sink), std::ref(graph));

            //t1.join();
            //t2.join();
            found = FoundMeet;//FoundMeet.load();
        }
        if (found) {
            meet_node = MeetnodeAtomic; //MeetnodeAtomic.load();
            std::cout << "Path found, meet node: " << meet_node << "\n";
            return true;
        } else {
            std::cout << "No path found!\n";
            return false;
        }   
    }

    void forward_search(std::vector<int>& parent, std::vector<int>& parent_rev, std::deque<int>& q_fwd, std::deque<int>& q_bwd, const std::vector<edges_matrix>& graph) {
        if (!q_fwd.empty()) {   //&&NOT FOUND
                int u = q_fwd.front();  
                q_fwd.pop_front();

                for (const auto& [to, edgeIdx] : FastNameMap[u]) {
                    const edges_matrix& a = graph[edgeIdx];
                    if (a.capacity - a.flow > err && parent[to] == -1) {    //&&NOT FOUND
                        parent[to] = u; //mtx lock
                        q_fwd.push_back(to);    //mtx lock
                                                //start new thread if pushed_back.??
                        if (parent_rev[to] != -1) {
                            q_fwd.clear();  //mtx lock
                            q_bwd.clear();  //mtx lock
                            MeetnodeAtomic = to; //MeetnodeAtomic.store(to);    
                            FoundMeet = true; //FoundMeet.store(true); 
                        }
                    }
                }
                for (int idx : reverse_adj[u]) {
                    int from = ReverseNameMapF[idx];
                    if (parent[from] == -1 && graph[idx].flow > err) {  //&&NOT FOUND
                        parent[from] = ReverseNameMapS[idx];    //mtx lock
                        q_fwd.push_back(from);  //mtx lock
                        if (parent_rev[from] != -1) {
                            q_fwd.clear();  //mtx lock
                            q_bwd.clear();  //mtx lock
                            MeetnodeAtomic = from; //MeetnodeAtomic.store(from);    
                            FoundMeet = true; //FoundMeet.store(true); 
                        }
                    }
                }
            }
    }

    void backward_search(std::vector<int>& parent, std::vector<int>& parent_rev, std::deque<int>& q_fwd, std::deque<int>& q_bwd, const int& sink, const std::vector<edges_matrix>& graph) {
        // Expand backward search
        if(!q_bwd.empty()) {
                int v = q_bwd.front();
                q_bwd.pop_front();

                for (int idx : reverse_adj[v]) {
                    int from = ReverseNameMapF[idx];
                    const edges_matrix& a = graph[idx];
                    if (a.capacity - a.flow > err && parent_rev[from] == -1) {  //&&NOT FOUND
                        parent_rev[from] = ReverseNameMapS[idx];    //mtx lock
                        q_bwd.push_back(from);  //mtx lock
                        if (parent[from] != -1) {
                            q_fwd.clear();  //mtx lock
                            q_bwd.clear();  //mtx lock
                            MeetnodeAtomic = from; //MeetnodeAtomic.store(from);    
                            FoundMeet = true; //FoundMeet.store(true); 
                        }
                    }
                }

                if(v != sink) {
                    for (const auto& [to, edgeIdx] : FastNameMap[v]) {
                        //if(to != sink){
                            if (graph[edgeIdx].flow > err && parent_rev[to] == -1) {  // Reverse in backward search //&&NOT FOUND
                                parent_rev[to] = v; //mtx lock
                                q_bwd.push_back(to);    //mtx lock
                                if (parent[to] != -1) {
                                    q_fwd.clear();  //mtx lock
                                    q_bwd.clear();  //mtx lock
                                    MeetnodeAtomic = to; //MeetnodeAtomic.store(to);    
                                    FoundMeet = true; //FoundMeet.store(true); 
                                }
                            }
                        //}
                    }
            }
        }
    }
*/
    void Edmonds_Karp() {
        const int source = sX;
        const int sink = tX;
        flow_val = 0;
        std::vector<int> parent(NameMap.size(), -1);
        parent[source] = source;

        //bool bfstool  = true;

        while (bfs(parent, source, sink, G_rk)) {

            /*clock_t z1 = clock();
            bfstool = bfs(parent, source, sink, G_rk);
            if(!bfstool) break;
            clock_t z2 = clock();
            totalBFS += (double)(z2-z1);*/

            bool found = false;
            double path_flow = std::numeric_limits<double>::max();

            // Find bottleneck
            for (int v = sink; v != source; v = parent[v]) {
                bool Reverse = false;
                int idx = NameMap[parent[v]][v];
                if (idx == 0) {
                    Reverse = true;
                    idx = NameMap[v][parent[v]];
                }
                edges_matrix x = G_rk[idx];

                if (!Reverse) {
                    if (x.capacity - x.flow > err) {
                        path_flow = std::min(path_flow, x.capacity - x.flow);
                        found = true;
                    }
                } else {
                    if (x.flow > err) {
                        path_flow = std::min(path_flow, x.flow);
                        found = true;
                    }
                }
            }

            // Update residual capacities
            for (int v = sink; v != source; v = parent[v]) {
                bool Reverse = false;
                int idx = NameMap[parent[v]][v];
                if (idx == 0) Reverse = true;

                if (!Reverse) {
                    if (G_rk[idx].capacity - G_rk[idx].flow > err) {
                        G_rk[idx].flow += path_flow;
                    }
                } else {
                    idx = NameMap[v][parent[v]];
                    G_rk[idx].flow -= path_flow;
                }
            }
        }
        // Compute total flow from I_a nodes to sink
        for (int i : I_a) {
            int Ikey = iX + i;
            flow_val += G_rk[NameMap[Ikey][tX]].flow;
        }
    }

    double totalBFS = 0;

    void Edmonds_Karp_Bidirectional() {
    const int source = sX;
    const int sink = tX;
    flow_val = 0;

    std::vector<int> parent(NameMap.size(), -1);
    std::vector<int> parent_rev(NameMap.size(), -1);
    int meet_node = -1;

    //bool bfstool = true;

    while (bidirectional_bfs(parent, parent_rev, meet_node, source, sink, G_rk)) {
        double path_flow = std::numeric_limits<double>::max();
        bool found = false;

        /*clock_t z1 = clock(); 
        bfstool = bidirectional_bfs(parent, parent_rev, meet_node, source, sink, G_rk);
        clock_t z2 = clock(); 
        totalBFS += (double)(z2-z1);
        if(!bfstool) break;*/

        // Reconstruct forward path: source -> meet_node
        for (int v = meet_node; v != source; v = parent[v]) {
            bool Reverse = false;
            int idx = NameMap[parent[v]][v];
            if (idx == 0) {
                Reverse = true;
                idx = NameMap[v][parent[v]];
            }
            const edges_matrix& x = G_rk[idx];

            if (!Reverse) {
                if (x.capacity - x.flow > err) {
                    path_flow = std::min(path_flow, x.capacity - x.flow);
                    found = true;
                }
            } else {
                if (x.flow > err) {
                    path_flow = std::min(path_flow, x.flow);
                    found = true;
                }
            }
        }
        
        // Reconstruct backward path: meet_node -> sink (via parent_rev)
        for (int v = meet_node; v != sink; v = parent_rev[v]) {
            int u = parent_rev[v];
            bool Reverse = false;
            int idx = NameMap[v][u];
            if (idx == 0) {
                Reverse = true;
                idx = NameMap[u][v];
            }
            const edges_matrix& x = G_rk[idx];

            if (!Reverse) {
                if (x.capacity - x.flow > err) {
                    path_flow = std::min(path_flow, x.capacity - x.flow);
                    found = true;
                }
            } else {
                if (x.flow > err) {
                    path_flow = std::min(path_flow, x.flow);
                    found = true;
                }
            }
        }
        

        if (!found) break;

        // Update backward path
        for (int v = meet_node; v != sink; v = parent_rev[v]) {
            int u = parent_rev[v];
            int idx = NameMap[v][u];
            if (G_rk[idx].capacity - G_rk[idx].flow > err) {
                G_rk[idx].flow += path_flow;
            } else {
                idx = NameMap[u][v];
                G_rk[idx].flow -= path_flow;
            }
        }

        // Update forward path
        for (int v = meet_node; v != source; v = parent[v]) {
            int u = parent[v];
            int idx = NameMap[u][v];
            if (G_rk[idx].capacity - G_rk[idx].flow > err) {
                G_rk[idx].flow += path_flow;
            } else {
                idx = NameMap[v][u];
                G_rk[idx].flow -= path_flow;
            }
        }
    }
    
    // Compute total flow from I_a nodes to sink
    for (int i : I_a) {
        int Ikey = iX + i;
        flow_val += G_rk[NameMap[Ikey][tX]].flow;
    }
}
    
    inline double length_sum_intervals(const std::vector<int>& I_list, const std::vector<int>& L_list) {
        double sum = 0.0;  // Initialize sum to 0
        for (int i : I_list) {  // Loop over each index in I_list
            sum += L_list[i];   // Add the corresponding value in L_list to sum
        }
        return sum;  // Return the total
    }

    double Get_M(std::vector<edges_matrix> X) {
            double M = 0;
            for (int j : jobs) {    //jobs is correct.
                int jobKey = jX + j;
                if (M == 0) {
                    M = X[NameMap[sX][jobKey]].capacity;
                }
                else {
                    M += X[NameMap[sX][jobKey]].capacity;
                }
            }
            return M;
    }

    int timeBase = 3600;
    std::vector<int> jobs;

    std::unordered_map<int, std::vector<int>> J;
    std::unordered_map<int, std::vector<int>> J_inverse;

    
    std::vector<int> intervals_start;
    std::vector<int> intervals_end;
    std::vector<int> I_a, I_p, I_crit_r;
    std::vector<int> len_i;
    std::vector<std::vector<int>> I_crit;

    std::vector<bool> subCrit_mask;
    std::vector<int> subCrit;

    std::vector<edges_matrix> G, G_r, G_rk, f_matrix;
    std::vector<std::map<int, int>> NameMap;
    std::vector<int> ReverseNameMapF, ReverseNameMapS;
    std::vector<std::vector<NameMapHelp>> FastNameMap;
    std::vector<std::vector<int>> reverse_adj;


    double EDMONDSKARPTIME = 0;
    int timestep;
    int rd, it;
    double flow_val = 0, flow_val_saved = 0, MaxDiff = 0;

    int sX = 0, jX = 1, tX = 0, iX = 0;

    double err = 0.0000001;
    bool selfterminate;
    double total_demand_r = 0;
    double objNormalized;
    bool AddingF = false;
};

#endif