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
#include <variant>
#include <cctype>
#include <algorithm>
#include <random>
#include <set>
#include <queue>

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

struct edges {
    std::string from;
    std::string to;
    double capacity;
    double flow;

    // Constructor
    edges() : from(""), to(""), capacity(0), flow(0) {}
    edges(std::string f, std::string t, double c, double fl) : from(f), to(t), capacity(c), flow(fl) {}
};

class Graph {
    public:
    std::vector<edges> digraph;

    Graph(int numnodes) : digraph(numnodes) {}

    void remove_empty(std::vector<edges>& digraph){
        for (int i =0; i < digraph.size(); ) {
            edges e = digraph.at(i);
            if (e.from == "" && e.to == "") {
                digraph.erase(digraph.begin() + i); //dont increase i, since we already move a slot.
            }
            else {
                i++;
            }
        }
    }

    //accept both an edge and simply the variables to make an edge.

    void add_flow(const edges& edge) {
        digraph.push_back(edge);
    }

    void add_flow(const std::string& from, const std::string& to, double capacity) {
        edges edge(from, to, capacity, 0);  // default flow is 0
        digraph.push_back(edge);
    }
    
    void print_graph(std::vector<edges> digraph) {
        for(int i =0; i < digraph.size(); i++) {
            edges e = digraph.at(i);
            std::cout << "From " << e.from << " to " << e.to
                          << " | cap: " << e.capacity << " | flow: " << e.flow << "\n";
        }
    }

    void find_J() {
        for (int i : I_a) {
            std::string i_key = "i" + std::to_string(i);
            std::vector<int> related_jobs;
        
            for (int j : jobs) {
                std::string j_key = "j" + std::to_string(j);
                const auto& job_list = J_inverse[j_key];
                if (std::find(job_list.begin(), job_list.end(), i) != job_list.end()) {
                    related_jobs.push_back(j);
                }
            }
        
            J[i_key] = related_jobs;
        }        
    }

    edges& GetEdge(std::string from, std::string to, std::vector<edges>& WhatGraph) {
        for (int i = 0; i < WhatGraph.size(); i++) {
            if(WhatGraph[i].from == from && WhatGraph[i].to == to) {
                return WhatGraph[i];
            }
        }
        std::cout << "from : " << from << " to : " << to << "\n";
        throw std::runtime_error("Edge not found");
    }

    void find_J_inverse() {
        for (int j = 0; j < jobs.size(); j++) {
            std::string key = "j" + std::to_string(j);
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
    
    void init_focs(InstanceData instance, int timeStep, int instancesize, bool randomize, int I_a_count) {

        remove_empty(digraph);

        timestep = timeStep;
        //init all

        average_power_W = sample_vector_C(instance.get_double_array("average_power_W"), instancesize, randomize);
        total_energy = sample_vector_C(instance.get_double_array("total_energy"), instancesize, randomize);
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
        find_J_inverse();

        //focs_instance_to_network
        for (int j : jobs) {
            //std::cout << j << "\n";
            std::string from = "s";
            std::string to = "j" + std::to_string(j);
            double capacity = jobs_demand[j];
            add_flow(from, to, capacity);  // D_s
        }
        
        for (int j : jobs) {
            std::string from = "j" + std::to_string(j);
            double cap_j = jobs_cap[j];
            std::string j_key = "j" + std::to_string(j);
        
            for (int i : J_inverse[j_key]) {
                std::string to = "i" + std::to_string(i);
                double capacity = cap_j * len_i[i] / timeBase;
                add_flow(from, to, capacity);  // D_0
            }
        }
        
        for (int i : I_a) {
            std::string from = "i" + std::to_string(i);
            std::string to = "t";
            double capacity = 0.0;  // If no value is given, use 0 or determine based on context
            add_flow(from, to, capacity);  // D_t
        }     
        digraph_r = digraph;
        total_demand_r = Get_M(digraph_r);  
        
        selfterminate = false;
        it = 0;
    }
    
    void reduce_network(std::vector<int> crit_r, std::vector<edges>& graph_rK) {
        for (int i : crit_r) {
            std::string Ikey = "i" + std::to_string(i);
            for (int j = 0; j < graph_rK.size(); j++) {
                if (graph_rK[j].from == Ikey || graph_rK[j].to == Ikey) {
                    graph_rK[j].from = "";
                    graph_rK[j].to = "";
                }
            }
        }
        remove_empty(graph_rK);
    }

    std::vector<edges> combineflows(std::vector<edges> start, std::vector<edges> adder) {
        for (int i = 0; i < start.size(); i++) {
            for (int j = 0; j < adder.size(); j++) {
                if (start[i].from == adder[j].from && start[i].to == adder[j].to) {
                    start[i].flow += adder[j].flow;
                    //std::cout << start[i].flow << " " << adder[j].flow << "\n";
                    break;
                }
            }
        }
        return start;
    }
    
    void partial_flow_func(std::vector<int> crit_r){ 
        partial_flow = digraph_rk;
        for (int i : I_a) {
            if (std::find(crit_r.begin(), crit_r.end(), i) == crit_r.end())  {  //is not in crit_r
                std::string Ikey = "i" + std::to_string(i);
                GetEdge(Ikey, "t", partial_flow).flow = 0;
                for (int j : J[Ikey]) {
                    std::string Jkey = "j" + std::to_string(j);
                    GetEdge(Jkey, Ikey, partial_flow).flow = 0;
                }
            }
        } 
        for (int j : jobs) {
            std::string Jkey = "j" + std::to_string(j);
            double sum_flow = 0.0;

            for (const edges& e : partial_flow) {
                if (e.from == Jkey) {
                    sum_flow += e.flow;
                }
            }
            GetEdge("s", Jkey, digraph_rk).flow = sum_flow;
        }
    }

    void update_network_capacities_g() {
        double demand = 0;
        if (it > 0) {
            demand = MaxDiff;
        }
        else {
            demand = Get_M(digraph_rk);
        }

        double demand_normalized = demand / length_sum_intervals(I_a, len_i);
        std::cout << "demand C++ = " << demand << "\n";
        std::cout << "length_sum_intervals(I_a, len_i) C++ = " << length_sum_intervals(I_a, len_i) << "\n";

        std::cout << "I_a ";
        for (int i = 0; i < I_a.size(); i++ ) {
            std::cout << I_a[i] << ", ";
        }
        std::cout << "\n";
        std::cout << "demand_normalized C++ = " << demand_normalized << "\n\n\n";

        for (int i : I_a) {
            std::string Ikey = "i" + std::to_string(i);
            GetEdge(Ikey, "t", digraph_rk).capacity = demand_normalized * len_i[i];
        }
    }

    void solve_focs() {
        while (!selfterminate) {
            if(it == 0) {
                digraph_rk = digraph_r;
            }
            update_network_capacities_g();
            
            Edmonds_Karp();
            MaxDiff = Get_M(digraph_rk) - flow_val; 

            if (total_demand_r-flow_val < err) {
                //end round

                //Check for subcrit. For some instances, there are still active intervals that only now don't reach the max anymore
                subCrit_mask.clear();
                std::string toKey = "t";
                for (int i : I_a) {
                    std::string fromKey = "i" + std::to_string(i);
                    std::cout << "fromKey: " << fromKey << "\n";
                    edges TheEdge = GetEdge(fromKey, toKey, digraph_rk);
                    /*std::cout << "Edge from " << TheEdge.from << " to " << TheEdge.to << " - "
              << "Capacity: " << TheEdge.capacity << ", Flow: " << TheEdge.flow
              << ", Err: " << err << std::endl;*/

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
                I_p = subCrit;
                std::cout << "I_p.size() = " << I_p.size() << " ////////////////////////////////////////////// \n";

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
                    f = combineflows(digraph_r, digraph_rk);
                }
                else {
                    I_crit_r = I_crit.back();
                    partial_flow_func(I_crit_r);
                    f = combineflows(partial_flow, digraph_rk);
                
                    reduce_network(I_crit_r, partial_flow);

                    I_a = I_p;                     // Copy the contents of I_p to I_a
                    std::sort(I_a.begin(), I_a.end());  // Sort I_a in ascending order
                    I_p.clear();

                    total_demand_r = Get_M(digraph_r);
                
                    rd++;
                    it = 0;
                }
            }


            else {
                
                subCrit_mask.clear();
                std::string toKey = "t";
                for (int i : I_a) {
                    std::string fromKey = "i" + std::to_string(i);
                    edges TheEdge = GetEdge(fromKey, toKey, digraph_rk);
                    bool isCritical = (TheEdge.capacity - TheEdge.flow > err);
                    subCrit_mask.push_back(isCritical);
                }

                subCrit.clear();
                for (size_t i = 0; i < I_a.size(); ++i) {
                    if (subCrit_mask[i]) {
                        subCrit.push_back(I_a[i]);
                    }
                }

                partial_flow_func(subCrit);
                reduce_network(subCrit, digraph_rk);

                total_demand_r = Get_M(digraph_rk);

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
        if (it <= 10) {
            objective();
        }
    }

    void Edmonds_Karp() {
        std::string source = "s";
        std::string sink = "t";
        flow_val = 0;
        std::unordered_map<std::string, std::string> parent;
    
        while (bfs(parent, source, sink, digraph_rk)) {
            // Find bottleneck capacity
            double path_flow = std::numeric_limits<double>::max();
            for (std::string v = sink; v != source; v = parent[v]) {
                std::string u = parent[v];
                
                // Check both forward and reverse edges for available capacity
                bool forward_found = false;
                bool reverse_found = false;
                
                for (edges& e : digraph_rk) {
                    if (e.from == u && e.to == v) {
                        // Forward edge: residual capacity is capacity - flow
                        if (e.capacity - e.flow > err) {
                            forward_found = true;
                            path_flow = std::min(path_flow, e.capacity - e.flow);
                        }
                    } else if (e.from == v && e.to == u) {
                        // Reverse edge: residual capacity is just the flow on the reverse edge
                        if (e.flow > err) {  // This assumes you are pushing flow through the reverse direction
                            reverse_found = true;
                            path_flow = std::min(path_flow, e.flow);
                        }
                    }
                }
            
                if (!forward_found && !reverse_found) {
                    std::cerr << "Error: No augmenting path found between " << u << " and " << v << std::endl;
                }
            }
    
            // Update residual capacities
            for (std::string v = sink; v != source; v = parent[v]) {
                std::string u = parent[v];
                // Forward edge
                bool forward_found = false;
                for (edges& e : digraph_rk) {
                    if (e.from == u && e.to == v) {
                        e.flow += path_flow;
                        forward_found = true;        
                        break;
                    }
                }
                // Reverse edge
                if (!forward_found) {
                    for (edges& e : digraph_rk) {
                        if (e.from == v && e.to == u) {
                            e.flow -= path_flow;
                            break;
                        }
                    }
                }

            }

            flow_val += path_flow;
        }
    }
    

    //Breadth First Search
    bool bfs(std::unordered_map<std::string, std::string>& parent, const std::string& source, const std::string& sink, const std::vector<edges>& graph) {
        std::unordered_map<std::string, bool> visited;
        std::queue<std::string> q;
        q.push(source);
        visited[source] = true;
        parent.clear();
    
        while (!q.empty()) {
            std::string u = q.front();
            q.pop();
            
            for (const edges& e : graph) {
                // Forward edge
                if (e.from == u && !visited[e.to] && e.capacity - e.flow > err) {
                    q.push(e.to);
                    parent[e.to] = u;
                    visited[e.to] = true;
                    if (e.to == sink) return true;
                }
                // Reverse edge
                else if (e.to == u && !visited[e.from] && e.flow > err) {
                    q.push(e.from);
                    parent[e.from] = u;
                    visited[e.from] = true;
                }
            }
            
            
        }
        return false;
    }

    void objective() {
        print_graph(f);
        int m = intervals_start.size();
        std::vector<double> p_i(m);
        for (int i = 0; i < m; i++) {
            std::string Ikey = "i" + std::to_string(i);
            double flow = GetEdge(Ikey, "t", f).flow;
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

    inline double length_sum_intervals(const std::vector<int>& I_list, const std::vector<int>& L_list) {
        double sum = 0.0;  // Initialize sum to 0
        for (int i : I_list) {  // Loop over each index in I_list
            sum += L_list[i];   // Add the corresponding value in L_list to sum
        }
        return sum;  // Return the total
    }

    double Get_M(std::vector<edges> X) {
            double M = 0;
            for (int j : jobs) {    //jobs is correct.
                std::string jobKey = "j" + std::to_string(j);
                if (M == 0) {
                    M = GetEdge("s", jobKey, X).capacity;
                }
                else {
                    M += GetEdge("s", jobKey, X).capacity;
                }
            }
            return M;
    }

    std::unordered_map<std::string, std::vector<int>> J;
    std::unordered_map<std::string, std::vector<int>> J_inverse;

    int timeBase = 3600;
    std::vector<int> jobs;
    std::vector<double> jobs_cap;
    std::vector<double> jobs_demand;
    std::vector<int> jobs_departure;
    std::vector<int> jobs_arrival;
    int n;

    std::vector<int> breakpoints;
    std::vector<int> intervals_start;
    std::vector<int> intervals_end;
    std::vector<int> I_a, I_p, I_crit_r;
    std::vector<int> len_i;
    std::vector<std::vector<int>> I_crit;

    std::vector<bool> subCrit_mask;
    std::vector<int> subCrit;

    std::vector<double> average_power_W;
    std::vector<double> t0_x;
    std::vector<double> t1_x;
    std::vector<double> total_energy;
    std::vector<double> total_energy_Wh;
    std::vector<double> maxPower;

    int timestep;
    int rd, it;
    std::vector<edges> digraph_r, digraph_rk, partial_flow, f;
    double flow_val = 0, MaxDiff = 0;

    double err = 0.0000001;
    bool MPCstopper = false;
    int MPCcondition = 0;
    bool selfterminate;
    double total_demand_r = 0;
    double objNormalized;
};

#endif