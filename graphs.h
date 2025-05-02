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

    void remove_empty(){
        for (int i =0; i < digraph.size(); ) {
            edges e = digraph.at(i);
            if (e.from.empty() && e.to.empty()) {
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
    
    void print_graph() {
        for(int i =0; i < digraph.size(); i++) {
            edges e = digraph.at(i);
            std::cout << "From " << e.from << " to " << e.to
                          << " | cap: " << e.capacity << "\n";
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

        timestep = timeStep;
        //init all
        for (int j = 0; j < I_a_count - 1; j++) {
            I_a[j] = j;
        }

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

        for (size_t i = 0; i < breakpoints.size() - 1; ++i) {
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
        
        selfterminate = false;
        digraph_r = digraph;
    }

    void update_network_capacities_g() {
        if (it > 0) {

        }
        else {

        }
    }
    
    void reduce_network() {
    
    }
    
    void length_sum_intervals() {
    
    } 
    
    void partial_flow(){ 
    
    }

    void solve_focs() {
        while (!selfterminate) {
            if(it == 0) {
                digraph_rk = digraph_r;
            }
            update_network_capacities_g();
            //flow_val = shortest_augmenting_path();
            std::cout << "flow_val = " << flow_val << "\n";
            selfterminate = true;
        }
    }
    
    double calculate_total_demand_r(){ 
        return 0;
    }

    private: 
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
    std::vector<int> I_a;
    std::vector<int> len_i;

    std::vector<double> average_power_W;
    std::vector<double> t0_x;
    std::vector<double> t1_x;
    std::vector<double> total_energy;
    std::vector<double> total_energy_Wh;
    std::vector<double> maxPower;

    int timestep;
    int rd, it;
    std::vector<edges> digraph_r;
    std::vector<edges> digraph_rk;
    double flow_val;

    double err = 0.0000001;
    bool MPCstopper = false;
    int MPCcondition = 0;
    bool selfterminate;
};

#endif