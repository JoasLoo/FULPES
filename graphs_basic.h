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
#include <unordered_set>
#include <functional>


inline std::vector<size_t> get_sample_indices(size_t dataSize, int N, bool randomSample) {
    std::vector<size_t> indices(dataSize);
    std::iota(indices.begin(), indices.end(), 0); // Fill with 0, 1, ..., n-1

    if (randomSample) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(indices.begin(), indices.end(), gen);
    }

    if (N < indices.size()) {
        indices.resize(N); // Keep only N samples
    }

    return indices;
}

inline std::vector<double> sample_by_indices(const std::vector<double>& input, const std::vector<size_t>& indices) {
    std::vector<double> result;
    result.reserve(indices.size());

    for (size_t idx : indices) {
        result.push_back(input[idx]);
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
    
    void init_focs(InstanceData instance, int timeStep, int instancesize, bool randomize) {

        remove_empty(digraph);

        timestep = timeStep;
        //init all

        // Get base vectors
        std::vector<double> average_power_base = instance.get_double_array("average_power_W");
        std::vector<double> total_energy_base = instance.get_double_array("total_energy_Wh");
        std::vector<double> max_power_base    = instance.get_double_array("maxPower");

        char t0_name[10], t1_name[10];
        snprintf(t0_name, sizeof(t0_name), "t0_%d", timestep);
        snprintf(t1_name, sizeof(t1_name), "t1_%d", timestep);
        std::vector<double> t0_base = instance.get_double_array(t0_name);
        std::vector<double> t1_base = instance.get_double_array(t1_name);

        // Shared indices for sampling
        std::vector<size_t> sample_indices = get_sample_indices(average_power_base.size(), instancesize, randomize);

        // Apply shared indices to all vectors
        average_power_W = sample_by_indices(average_power_base, sample_indices);
        total_energy_Wh = sample_by_indices(total_energy_base, sample_indices);
        maxPower        = sample_by_indices(max_power_base, sample_indices);
        t0_x            = sample_by_indices(t0_base, sample_indices);
        t1_x            = sample_by_indices(t1_base, sample_indices);
        
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
            double remove_from_s_to_jX = 0;
            for (int j = 0; j < graph_rK.size(); j++) {
                if (graph_rK[j].from == Ikey || graph_rK[j].to == Ikey) {
                    if(graph_rK[j].to == Ikey) {
                        remove_from_s_to_jX = graph_rK[j].flow;
                        GetEdge("s", graph_rK[j].from, digraph_rk).capacity -= remove_from_s_to_jX;
                        GetEdge("s", graph_rK[j].from, digraph_rk).flow -= remove_from_s_to_jX;
                    }
                    if(AddingF) {  // if F is initiated in another round, add the flows to F.
                        GetEdge(graph_rK[j].from, graph_rK[j].to, f).flow += graph_rK[j].flow;
                    }
                    graph_rK[j].from = "";
                    graph_rK[j].to = "";
                }
            }
        }
        remove_empty(graph_rK);
    }



    void Add_to_f_end(std::vector<edges>& graph_rK) {
        for (int i = 0; i < graph_rK.size(); i++) {
            GetEdge(graph_rK[i].from, graph_rK[i].to, f).flow += graph_rK[i].flow;
        }
    }

    std::vector<edges> combineflows(std::vector<edges> start, std::vector<edges> adder) {
        if(start.empty()) {
            return adder;
        }
        for (int i = 0; i < start.size(); i++) {
            for (int j = 0; j < adder.size(); j++) {
                if (start[i].from == adder[j].from && start[i].to == adder[j].to) {
                    start[i].flow += adder[j].flow;
                    break;
                }
            }
        }
        return start;
    }
    
    void partial_flow_func(std::vector<int> crit_r){ 
        partial_flow = digraph_r;
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
            GetEdge("s", Jkey, digraph_r).flow = sum_flow;
        }
    }

    void reset_flows(std::vector<edges>& GRAPH) {
        for (int i = 0; i < GRAPH.size(); i++) {
            GRAPH[i].flow = 0;
        } 
    }

    void reset_caps(std::vector<edges>& GRAPH) {
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
                std::string Ikey = "i" + std::to_string(i);
                GetEdge(Ikey, "t", digraph_rk).capacity += demand_normalized * len_i[i];
            }
        }
        else {
            if (rd == 0) {
                demand = Get_M(digraph_rk);
            }
            else {
                demand = Get_M(digraph_rk)-flow_val_saved;
            }
            demand_normalized = demand / length_sum_intervals(I_a, len_i);
            for (int i : I_a) {
                std::string Ikey = "i" + std::to_string(i);
                GetEdge(Ikey, "t", digraph_rk).capacity = demand_normalized * len_i[i];
            }
        }
    }

    void solve_focs() {
        f = digraph;
        reset_caps(f);
        reset_flows(f);
        while (!selfterminate) {
            if(it == 0) {
                digraph_rk = digraph_r;
            }
            update_network_capacities_g();
            //Edmonds_Karp();
            Max_flow_solver();
            MaxDiff = Get_M(digraph_rk) - flow_val; 
            
            if (total_demand_r-flow_val < err) {
                //end round

                //Check for subcrit. For some instances, there are still active intervals that only now don't reach the max anymore
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
                    AddingF = true; //we need this otherwise F gets added also when subcrit. but we only want when crit.
                    Add_to_f_end(digraph_rk);
                    AddingF = false;    //still setting to false so i dont have to do that in init_Focs anymore <3
                }
                else {
                    I_crit_r = I_crit.back();
                    AddingF = true;
                    reduce_network(I_crit_r, digraph_rk);
                    AddingF = false;
                    I_a = I_p;                     // Copy the contents of I_p to I_a
                    std::sort(I_a.begin(), I_a.end());  // Sort I_a in ascending order
                    I_p.clear();

                    total_demand_r = Get_M(digraph_r)-flow_val_saved;

                    flow_val_saved += flow_val;
                
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

                //partial_flow_func(subCrit);
                reduce_network(subCrit, digraph_rk);

                total_demand_r = Get_M(digraph_rk)-flow_val_saved;

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
    
    void Max_flow_solver() {
        Edmonds_Karp();
        // After redistributing, run Edmonds-Karp again to fill freed space
        for (int i : I_a) {
            std::string ikey = "i" + std::to_string(i);
            double cap_to_t = GetEdge(ikey, "t", digraph_rk).capacity;
            double flow_to_t = GetEdge(ikey, "t", digraph_rk).flow;
        
            if (flow_to_t + 0.0001 < cap_to_t) {
                reset_flows(digraph_rk);   // Reset all flows to 0
                Edmonds_Karp();            // Recompute full max-flow from scratch
                break;                     // Exit loop; re-evaluate from start after EK
            }
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
        }
        double temp;
        for (int i : I_a) {
            std::string Ikey = "i" + std::to_string(i);
            flow_val += GetEdge(Ikey, "t", digraph_rk).flow;
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

    bool redistribute_flow(const std::string& target_node, double needed_capacity) {
        std::unordered_map<std::string, std::string> parent;
        std::unordered_set<std::string> visited;
    
        std::function<bool(std::string)> dfs = [&](std::string u) {
            visited.insert(u);
            if (u == "s") return true;
    
            for (edges& e : digraph_rk) {
                // Traverse reverse edges (i.e., flow > 0)
                if (e.to == u && e.flow > err && visited.find(e.from) == visited.end()) {
                    parent[e.from] = u;
                    if (dfs(e.from)) return true;
                }
            }
            return false;
        };
    
        if (!dfs(target_node)) return false;
    
        // Find bottleneck along reverse path
        double path_flow = needed_capacity;
        for (std::string v = "s"; v != target_node; v = parent[v]) {
            std::string u = parent[v];
            for (edges& e : digraph_rk) {
                if (e.to == u && e.from == v) {
                    path_flow = std::min(path_flow, e.flow);
                }
            }
        }
    
        // Apply flow reduction along reverse path
        for (std::string v = "s"; v != target_node; v = parent[v]) {
            std::string u = parent[v];
            for (edges& e : digraph_rk) {
                if (e.to == u && e.from == v) {
                    e.flow -= path_flow;
                }
            }
        }
    
        return true;
    }    

    void objective() {
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
    double flow_val = 0, flow_val_saved = 0, MaxDiff = 0;

    double err = 0.0000001;
    bool MPCstopper = false;
    int MPCcondition = 0;
    bool selfterminate;
    double total_demand_r = 0;
    double objNormalized;
    bool AddingF = false;
};

#endif