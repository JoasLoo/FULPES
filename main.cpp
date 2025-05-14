#include "graphs.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace std;


const int instanceSize = 200; //number of EVs/jobs in instance
int timeStep = 900; //quarterly granularity

bool randomSample = false;

InstanceData opendata_toC(const string& filename);

int main() {
    FILE *data; 
    data = fopen("Data/DEMSdata_FOCS_v1.csv","r"); //open data file in C DEMSdata_FOCS_v1.csv ev_session_data_OR.csv
    if(data == NULL) { //Check if data file opened succesfully
        printf("Cant open data.csv");
        exit(0);
    }
    
    clock_t t1 = clock();
    InstanceData instance1 = opendata_toC("Data/DEMSdata_FOCS_v1.csv");
    clock_t q1 = clock();

    Graph g;
    clock_t q2 = clock();
    
    g.init_focs(instance1, timeStep, instanceSize, randomSample, 3600);
    clock_t q3 = clock(); 
    g.solve_focs();
    
    clock_t qx = clock();


    printf("X) LOADING C DATA            %.5f seconds\n", (double)(q1 - t1) / CLOCKS_PER_SEC);
    printf("X) INIT C GRAPH            %.5f seconds\n", (double)(q2 - q1) / CLOCKS_PER_SEC);
    printf("X) INIT C FOCS            %.5f seconds\n", (double)(q3 - q2) / CLOCKS_PER_SEC);
    printf("X) SOLVE C FOCS            %.5f seconds\n", (double)(qx - q3) / CLOCKS_PER_SEC);

    fclose(data);
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
            instance1.add_column(headers[i], move(col_d));
        }

        instance1.add_column(headers[i], move(datainCols[i]));
        //cout << headers[i] << " " << datainCols[i][0] << "\n";
    }
    
    file.close();
    return instance1;
}
