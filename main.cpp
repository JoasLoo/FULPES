#include "graphs.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
using namespace std;
using namespace std::chrono;

const int instanceSize = 2000; //number of EVs/jobs in instance
int timeStep = 900; //quarterly granularity

bool randomSample = false;

vector<int> FOCS_breakpoints;

InstanceData opendata_toC(const string& filename);


int main() {

    long long total_load = 0;
    long long total_graph = 0;
    long long total_init = 0;
    long long total_solve = 0;
    
    int repetitions = 100;


    auto t1 = chrono::high_resolution_clock::now();

    InstanceData instance1 = opendata_toC("Data/DEMSdata_FOCS_v1.csv"); //open data file in C++ DEMSdata_FOCS_v1.csv ev_session_data_OR.csv
        //ev_session_data_OR.csv breaks if randomsample, or instancesize = 399 / 200

    for (int i = 0; i < repetitions; ++i) {
    
    auto q1 = chrono::high_resolution_clock::now();
    Graph g;
    //Graph g(1);
    //g.remove_empty(g.digraph);
    auto q2 = chrono::high_resolution_clock::now();
    
    g.init_focs(instance1, timeStep, instanceSize, randomSample); //InstanceData instance, int timeStep, int instancesize, bool randomize, int I_a_count
    auto q3 = chrono::high_resolution_clock::now();
    //g.print_graph();   
    g.solve_focs();
    
    auto qx = chrono::high_resolution_clock::now();

    g.objective();

    if (i==0) {
        total_load  += duration_cast<microseconds>(q1 - t1).count();
    }
    
    total_graph += duration_cast<microseconds>(q2 - q1).count();
    total_init  += duration_cast<microseconds>(q3 - q2).count();
    total_solve += duration_cast<microseconds>(qx - q3).count();

    }
    cout << fixed;
    //cout << "Average Timings (over " << repetitions << " runs):\n";
    //cout << "LOADING C DATA       " << total_load  << " micro-s\n";
    //cout << "INIT C GRAPH         " << total_graph / repetitions << " micro-s\n";
    //cout << "INIT C FOCS          " << total_init / repetitions  << " micro-s\n";
    cout << "SOLVE C FOCS         " << total_solve / repetitions << " micro-s\n";



    return 0;
}

bool starts_with(const std::string& str, const std::string& prefix) {
    return str.size() >= prefix.size() &&
           std::equal(prefix.begin(), prefix.end(), str.begin());
}

InstanceData opendata_toC(const string& filename) {   //copied from: https://medium.com/@ryan_forrester_/reading-csv-files-in-c-how-to-guide-35030eb378ad

    InstanceData instance1;
    vector<vector<string>> data;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return instance1;
    }

    string line;
    while (getline(file, line)) {
        vector<string> row;
        stringstream ss(line);
        string cell;

        while (getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        data.push_back(row);
    }

    vector<vector<string>> datainCols;
    vector<string> headers;

    for (int i = 0; i < data[0].size(); i++) {
        vector<string> col;
        for (int q = 0; q < data.size(); q++) {
            if (q == 0) {
                headers.push_back(data[q][i]);
            }
            else {
                col.push_back(data[q][i]);
            }
            
        }
        datainCols.push_back(col);
    }

    if (headers[0] == "") {
        headers[0] = "Unnamed ";
    }

    for (int i = 0; i < data[0].size(); i++) { 
        vector<string>& col = datainCols[i];
        if(starts_with(headers[i], "t0_") || starts_with(headers[i], "t1_") ||
        headers[i] == "total_energy" || headers[i] == "total_energy_Wh" ||
        headers[i] == "average_power_W" || headers[i] == "maxPower") {
            vector<double> col_d;
            col_d.reserve(col.size());
            for (const auto& val : col) {
                try {
                    col_d.push_back(std::stod(val));
                } catch (const std::invalid_argument& e) {
                    col_d.push_back(0.0); // Or handle invalid entries differently
                }
            }
            instance1.add_column(headers[i], std::move(col_d));
        }

        instance1.add_column(headers[i], std::move(datainCols[i]));
        //cout << headers[i] << " " << datainCols[i][0] << "\n";
    }
    
    file.close();
    return instance1;
}
