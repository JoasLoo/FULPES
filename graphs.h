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
#include <variant>
#include <cctype>
#include <algorithm>
#include <random>
#include <set>
#include <queue>
#include <unordered_set>
#include <functional>
#include <iomanip>  // For std::setw and std::setprecision
#include <tuple>

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

class Graph {
    public:

    void print_graph(const std::vector<double>& digraph_cap, const std::vector<double>& digraph_flow, std::unordered_map<int, std::pair<int, int>> ReverseNameMap) {
        std::cout << "-------------------------------------------- \n";
    for (size_t i = 0; i < digraph_cap.size(); ++i) {
        double flow = digraph_flow[i];
        double capacity = digraph_cap[i];

        auto it = ReverseNameMap.find(i);
        if (it != ReverseNameMap.end()) {
            const int& from = it->second.first;
            const int& to = it->second.second;
            std::cout << "From: " << from << " To: " << to
                      << " | Capacity: " << capacity
                      << " | Flow: " << flow << "\n";
        }
    }
}


    void find_J(std::unordered_map<int, std::vector<int>>& J, std::unordered_map<int, std::vector<int>>& J_inverse, std::vector<int> I_a) {
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

    void find_J_inverse(std::unordered_map<int, std::vector<int>>& J_inverse, std::vector<double> t0_x, std::vector<double> t1_x, std::vector<int> intervals_end, std::vector<int> intervals_start) {
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
        //          G_r_capacity, G_r_flow, G_rk_capacity, G_rk_flow, f, I_a, NameMap, ReverseNameMap, reverse_adj 
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int>, std::unordered_map<int, std::map<int, int>>, std::unordered_map<int, std::pair<int, int>>, std::unordered_map<int, std::vector<int>>> init_focs(InstanceData instance, int timeStep, int instancesize, bool randomize, int timeBase) {

        std::vector<double> jobs_cap;
        std::vector<double> jobs_demand;
        std::vector<int> jobs_departure;
        std::vector<int> jobs_arrival;

        std::vector<int> breakpoints;
        std::vector<int> intervals_end;

        std::vector<double> average_power_W;
        std::vector<double> t0_x;
        std::vector<double> t1_x;
        std::vector<double> total_energy;
        std::vector<double> total_energy_Wh;
        std::vector<double> maxPower;

        std::unordered_map<int, std::vector<int>> J;
        std::unordered_map<int, std::vector<int>> J_inverse;

        std::vector<int> I_a;
        std::vector<double> G, G_r_capacity, G_r_flow, G_rk_capacity, G_rk_flow, f;
        std::unordered_map<int, std::map<int, int>> NameMap;
        std::unordered_map<int, std::pair<int, int>> ReverseNameMap;
        std::unordered_map<int, std::vector<int>> reverse_adj;

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

        find_J(J, J_inverse, I_a);
        find_J_inverse(J_inverse, t0_x, t1_x, intervals_end, intervals_start);

        //HERE INSTANTIATION OF MATRIX

        iX = jobs.size() + 1;
        tX = iX + I_a.size();
        //First add all names to NameMap and ReverseNameMap
        int edgeIndex = 1;      //Starting at 1 so if it defaults to 0 we know it couldnt find the edge.

        // Add all edge names to NameMap and ReverseNameMap
        for (int j : jobs) {
            int from = sX;
            int to = jX + j;
            NameMap[from][to] = edgeIndex;
            ReverseNameMap[edgeIndex].first = from;
            ReverseNameMap[edgeIndex].second = to;
            edgeIndex++;
        }

        for (int j : jobs) {
            int from = jX + j;
            double cap_j = jobs_cap[j];
            
            for (int i : J_inverse[from]) {
                int to = iX + i;
                NameMap[from][to] = edgeIndex;
                ReverseNameMap[edgeIndex].first = from;
            ReverseNameMap[edgeIndex].second = to;
                edgeIndex++;
            }
        }

        for (int i : I_a) {
            int from = iX + i;
            int to = tX;
            NameMap[from][to] = edgeIndex;
            ReverseNameMap[edgeIndex].first = from;
            ReverseNameMap[edgeIndex].second = to;
            edgeIndex++;
        }

        //dedicate enough space to all matrices.
        G = std::vector<double>(edgeIndex);
        G_rk_capacity = G;
        G_rk_flow = G;
        G_r_flow = G;
        f = G;
        //now we can access each matrix via: G[NameMap[sX]["j0"]].capacity = 1;


        //focs_instance_to_network
        for (int j : jobs) {
            int from = sX;
            int to = jX + j;
            double capacity = jobs_demand[j];
            G[NameMap[from][to]] = capacity;
        }
        
        for (int j : jobs) {
            int from = jX + j;
            double cap_j = jobs_cap[j];
        
            for (int i : J_inverse[from]) {
                int to = iX + i;
                double capacity = cap_j * len_i[i] / timeBase;
                G[NameMap[from][to]] = capacity;
            }
        }
        
        for (int i : I_a) {
            int from = iX + i;
            int to = tX;
            double capacity = 0.0;  // If no value is given, use 0 or determine based on context
            G[NameMap[from][to]] = capacity;
        }     

        for (const auto& [from, innerMap] : NameMap) {
            for (const auto& [to, edgeIdx] : innerMap) {
                reverse_adj[to].push_back(edgeIdx);  // Edge ends at `to`, so it’s a reverse edge for `to`
            }
        }

        G_r_capacity = G;

        total_demand_r = Get_M(G_r_capacity, NameMap);  
        selfterminate = false;
        it = 0;

        return std::make_tuple(G_r_capacity, G_r_flow, G_rk_capacity, G_rk_flow, f, I_a, NameMap, ReverseNameMap, reverse_adj);
    }

    void reduce_network(std::vector<int> crit_r, std::vector<double>& graph_rK_flow, std::vector<double>& graph_rK_capacity, bool AddingF, std::unordered_map<int, std::pair<int, int>> ReverseNameMap, std::unordered_map<int, std::vector<int>> reverse_adj, std::vector<double> f, std::unordered_map<int, std::map<int, int>> NameMap) {
        for (int i : crit_r) {
            int Ikey = iX + i;
            double remove_from_s_to_jX = 0;

            for (int j : reverse_adj[Ikey]) {   //if goes towards iX
                remove_from_s_to_jX = graph_rK_flow[j];
                graph_rK_capacity[NameMap[sX][ReverseNameMap[j].first]] -= remove_from_s_to_jX;
                graph_rK_flow[NameMap[sX][ReverseNameMap[j].first]] -= remove_from_s_to_jX;;
                if(AddingF) {  // if F is initiated in another round, add the flows to F.
                    f[j] += graph_rK_flow[j];
                }
            }

            for (const auto& [to, j] : NameMap[Ikey]) { //for goes from iX
                if(AddingF) {  // if F is initiated in another round, add the flows to F.
                    f[j] += graph_rK_flow[j];
                }
                graph_rK_capacity[j] = 0;
                graph_rK_flow[j] = 0;
            }
        }
    }

    void Add_to_f_end(std::vector<double>& graph_rK, std::vector<double>& f) {
        for (int i = 0; i < graph_rK.size(); i++) {
            f[i] += graph_rK[i];
        }
    }

    void reset(std::vector<double>& GRAPH) {
        for (int i = 0; i < GRAPH.size(); i++) {
            GRAPH[i] = 0;
        }
    }

    void update_network_capacities_g(std::vector<double>& G_rk, std::vector<int> I_a, std::unordered_map<int, std::map<int, int>> NameMap) {
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
                G_rk[NameMap[Ikey][tX]] += demand_normalized * len_i[i];
            }
        }
        else {
            if (rd == 0) {
                demand = Get_M(G_rk, NameMap);
            }
            else {
                demand = Get_M(G_rk, NameMap)-flow_val_saved;
            }
            demand_normalized = demand / length_sum_intervals(I_a, len_i);
            for (int i : I_a) {
                int Ikey = iX + i;
                G_rk[NameMap[Ikey][tX]] = demand_normalized * len_i[i];
            }

            //std::cout << "Demand: " << demand << " length_sum_intervals(I_a, len_i) " << length_sum_intervals(I_a, len_i) << "\n";
        }
    }


    void solve_focs(std::vector<double> G_r_capacity, std::vector<double> G_r_flow, std::vector<double> G_rk_capacity, std::vector<double> G_rk_flow, std::vector<double> f, std::vector<int> I_a, std::unordered_map<int, std::map<int, int>> NameMap, std::unordered_map<int, std::pair<int, int>> ReverseNameMap, std::unordered_map<int, std::vector<int>> reverse_adj) {
        reset(f);
        while (!selfterminate) {
            if(it == 0) {
                G_rk_capacity = G_r_capacity;
                G_rk_flow = G_r_flow;
            }
            update_network_capacities_g(G_rk_capacity, I_a, NameMap);
            //Edmonds_Karp();
            Max_flow_solver(G_rk_flow, G_rk_capacity, ReverseNameMap, NameMap, I_a, reverse_adj);
            MaxDiff = Get_M(G_rk_capacity, NameMap) - flow_val; 
            
            if (total_demand_r-flow_val < err) {
                //end round

                //Check for subcrit. For some instances, there are still active intervals that only now don't reach the max anymore
                subCrit_mask.clear();
                int toKey = tX;
                for (int i : I_a) {
                    int fromKey = iX + i;
                    double Capacity = G_rk_capacity[NameMap[fromKey][toKey]];
                    double Flow = G_rk_flow[NameMap[fromKey][toKey]];

                    bool isCritical = (Capacity - Flow > err);
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
                    Add_to_f_end(G_rk_flow, f);
                }
                else {
                    I_crit_r = I_crit.back();
                    AddingF = true;
                    reduce_network(I_crit_r, G_rk_flow, G_rk_capacity, AddingF, ReverseNameMap, reverse_adj, f, NameMap);
                    AddingF = false;
                    I_a = I_p;                     // Copy the contents of I_p to I_a
                    std::sort(I_a.begin(), I_a.end());  // Sort I_a in ascending order
                    I_p.clear();

                    total_demand_r = Get_M(G_r_capacity, NameMap)-flow_val_saved;

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
                    double Capacity = G_rk_capacity[NameMap[fromKey][toKey]];
                    double Flow = G_rk_flow[NameMap[fromKey][toKey]];

                    bool isCritical = (Capacity - Flow > err);
                    subCrit_mask.push_back(isCritical);
                }

                subCrit.clear();
                for (size_t i = 0; i < I_a.size(); ++i) {
                    if (subCrit_mask[i]) {
                        subCrit.push_back(I_a[i]);
                    }
                }

                reduce_network(subCrit, G_rk_flow, G_rk_capacity, AddingF, ReverseNameMap, reverse_adj, f, NameMap);

                total_demand_r = Get_M(G_rk_capacity, NameMap)-flow_val_saved;

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
            objective(3600, NameMap, f);
            printf("X) TOTAL EDMONDS KARP               %.5f seconds\n", EDMONDSKARPTIME / CLOCKS_PER_SEC); 
            printf("X) TOTAL BFS                        %.5f seconds\n", totalBFS / CLOCKS_PER_SEC); 
            printf("X) REST EDMONDS KARP                %.5f seconds\n", restEdmondsKarp / CLOCKS_PER_SEC); 
        }
    }
    
    void Max_flow_solver(std::vector<double> G_rk_flow, std::vector<double> G_rk_capacity, std::unordered_map<int, std::pair<int, int>> ReverseNameMap, std::unordered_map<int, std::map<int, int>> NameMap, std::vector<int> I_a, std::unordered_map<int, std::vector<int>> reverse_adj) {    //G_rk_flow, G_rk_capacity, ReverseNameMap, NameMap, I_a, reverse_adj
        clock_t z1 = clock();
        reset(G_rk_flow);   // Reset all flows to 0 for a DFS effect
        Edmonds_Karp(G_rk_flow, G_rk_capacity, ReverseNameMap, NameMap, I_a, reverse_adj);
        clock_t z2 = clock();
        EDMONDSKARPTIME += (double)(z2-z1);
    }
    double totalBFS = 0;
    double restEdmondsKarp = 0;

    void Edmonds_Karp(std::vector<double> G_rk_flow, std::vector<double> G_rk_capacity, std::unordered_map<int, std::pair<int, int>> ReverseNameMap, std::unordered_map<int, std::map<int, int>> NameMap, std::vector<int> I_a, std::unordered_map<int, std::vector<int>> reverse_adj) {
        int source = sX;
        int sink = tX;
        flow_val = 0;
        std::unordered_map<int, int> parent;

        int keepingcount = 0;

        bool bfstool = true;
        double TOTALBFS = 0;

        while (bfstool) {
            keepingcount++;
            bool found = false;
            double path_flow = std::numeric_limits<double>::max();
            clock_t z1 = clock();
            bfstool = bfs(parent, source, sink, G_rk_flow, G_rk_capacity, ReverseNameMap, NameMap, reverse_adj);
            if (!bfstool) break;
            clock_t z2 = clock();

            // Find bottleneck
            for (int v = sink; v != source; v = parent[v]) {
                bool Reverse = false;
                int idx = NameMap[parent[v]][v];
                if (idx == 0) {
                    Reverse = true;
                    idx = NameMap[v][parent[v]];
                }
                double x = G_rk_capacity[idx];
                double y = G_rk_flow[idx];

                if (!Reverse) {
                    if (x - y > err) {
                        path_flow = std::min(path_flow, x - y);
                        found = true;
                    }
                } else {
                    if (y > err) {
                        path_flow = std::min(path_flow, y);
                        found = true;
                    }
                }

                if (!found) {
                    std::cerr << "Error: No valid path from " << parent[v] << " to " << v << std::endl;
                }
            }

            // Update residual capacities
            for (int v = sink; v != source; v = parent[v]) {
                bool Reverse = false;
                int idx = NameMap[parent[v]][v];
                if (idx == 0) Reverse = true;

                if (!Reverse) {
                    if (G_rk_capacity[idx] - G_rk_flow[idx] > err) {
                        G_rk_flow[idx] += path_flow;
                    }
                } else {
                    idx = NameMap[v][parent[v]];
                    G_rk_flow[idx] -= path_flow;
                }
            }
            clock_t z3 = clock();
            totalBFS += (double)(z2-z1);
            restEdmondsKarp += (double)(z3-z2);
        }
        // Compute total flow from I_a nodes to sink
        for (int i : I_a) {
            int Ikey = iX + i;
            flow_val += G_rk_flow[NameMap[Ikey][tX]];
        }
    }

    //Breadth First Search
    bool bfs(std::unordered_map<int, int>& parent, const int& source, const int& sink, const std::vector<double> graph_flow, const std::vector<double> graph_capacity, std::unordered_map<int, std::pair<int, int>> ReverseNameMap, std::unordered_map<int, std::map<int, int>> NameMap, std::unordered_map<int, std::vector<int>> reverse_adj) {
        std::unordered_map<int, bool> visited;
        std::queue<int> q;
        q.push(source);
        visited[source] = true;
        parent.clear();

        int N = NameMap.size();
    
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            bool continue_outer = false;

            for (const auto& [to, edgeIdx] : NameMap[u]) {
                // `to` is the destination node name (e.g., "j0")
                // `edgeIdx` is the unique index for edge (sX -> to)

                if (!visited[to] && graph_capacity[edgeIdx] - graph_flow[edgeIdx] > err) {
                    q.push(to);
                    parent[to] = u;
                    visited[to] = true;
                    if (to == sink) return true;
                    continue_outer = true;
                }
            }
            if (continue_outer) continue;
            
            for (int idx : reverse_adj[u]) {
                const int& from = ReverseNameMap[idx].first;
                const int& to = ReverseNameMap[idx].second;
                if (to == u && !visited[from] && graph_flow[idx] > err) {
                    q.push(from);
                    parent[from] = to;
                    visited[from] = true;
                    continue_outer = true;
                }
            }
        }
        return false;
    }

    void objective(int timeBase, std::unordered_map<int, std::map<int, int>> NameMap, std::vector<double> f) {
        int m = intervals_start.size();
        std::vector<double> p_i(m);
        for (int i = 0; i < m; i++) {
            int Ikey = iX + i;
            double flow = f[NameMap[Ikey][tX]];
            p_i[i] = (flow / len_i[i]) * timeBase;
        }
        std::vector<double> powerSquare(m);
        for (int i = 0; i < m; ++i) {
            powerSquare[i] = (p_i[i] * p_i[i]) * (len_i[i] / timestep);
        }

        objNormalized = std::accumulate(powerSquare.begin(), powerSquare.end(), 0.0);
        std::cout << "\nObjective value C++ = " << objNormalized << "\n";
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

    double Get_M(std::vector<double> X, std::unordered_map<int, std::map<int, int>> NameMap) {
            double M = 0;
            for (int j : jobs) {    //jobs is correct.
                int jobKey = jX + j;
                if (M == 0) {
                    M = X[NameMap[sX][jobKey]];
                }
                else {
                    M += X[NameMap[sX][jobKey]];
                }
            }
            return M;
    }

    std::vector<int> jobs;
    int n;

    
    std::vector<int> intervals_start;

    std::vector<int> I_p, I_crit_r;
    std::vector<int> len_i;
    std::vector<std::vector<int>> I_crit;

    std::vector<bool> subCrit_mask;
    std::vector<int> subCrit;

    double EDMONDSKARPTIME = 0;
    int timestep;
    int rd, it;
    double flow_val = 0, flow_val_saved = 0, MaxDiff = 0;

    int sX = 0, jX = 1, tX = 0, iX = 0;

    double err = 0.0000001;
    bool MPCstopper = false;
    int MPCcondition = 0;
    bool selfterminate;
    double total_demand_r = 0;
    double objNormalized;
    bool AddingF = false;
};

#endif