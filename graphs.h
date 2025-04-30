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


class InstanceData {
public:
    std::unordered_map<std::string, std::vector<std::string>> data;

    // Access like: data["average_power_W"][j]
    std::string get(const std::string& key, size_t index) const {
        return data.at(key).at(index);
    }

    size_t size() const {
        if (data.empty()) return 0;
        return data.begin()->second.size();
    }

    void add_column(const std::string& name, const std::vector<std::string>& column) {
        data[name] = column;
    }
};




struct edges {
    std::string from;
    std::string to;
    double capacity;

    // Constructor
    edges() : from(""), to(""), capacity(0) {}
    edges(std::string f, std::string t, double c) : from(f), to(t), capacity(c) {}
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

    void add_flow(edges edge) {
        digraph.push_back(edge);
    }
    
    void print_graph() {
        for(int i =0; i < digraph.size(); i++) {
            edges e = digraph.at(i);
            std::cout << "From " << e.from << " to " << e.to
                          << " | cap: " << e.capacity << "\n";
        }
    }
    
    void init_focs() {
    
    }
    
    void reduce_network() {
    
    }
    
    void length_sum_intervals() {
    
    } 
    
    void partial_flow(){ 
    
    }
    
    double calculate_total_demand_r(){ 
        return 0;
    }
};

#endif