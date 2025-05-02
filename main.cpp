#include "link.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace std;


const int instanceSize = 70; //number of EVs/jobs in instance
int timeStep = 900; //quarterly granularity

//PyObject *maxFlowAlg = shortest_augmenting_path;  // #alternatively use e.g., edmonds_karp, preflow_push, or dinitz
bool randomSample = false;

PyObject *empty_flow_func = NULL;
PyObject *calculate_total_demand_r_func = NULL;

vector<int> FOCS_breakpoints;


void Fail_EXIT(const char *msg)
{
    fprintf(stderr, "%s\n", msg);
    PyErr_Print();
    Py_Finalize();
    exit(EXIT_FAILURE);
}

struct DataHelp {
    PyObject *Cols;
    PyObject *Names;
};


struct DataHelp getDataC(FILE *data);
InstanceData opendata_toC(const string& filename);
int extract_unique_sorted_times(struct DataHelp data);


int main() {
    FILE *data; 
    data = fopen("Data/DEMSdata_FOCS_v1.csv","r"); //open data file in C DEMSdata_FOCS_v1.csv ev_session_data_OR.csv
    if(data == NULL) { //Check if data file opened succesfully
        printf("Cant open data.csv");
        exit(0);
    }
    //Start Python environment
    Py_Initialize();
    
    //Import FOCS.py
    PyObject *pName = PyUnicode_DecodeFSDefault("FOCS");
    PyObject *FOCS = PyImport_Import(pName);
    Py_DECREF(pName);
    
    //Import script.py
    PyObject *pName1 = PyUnicode_DecodeFSDefault("script");
    PyObject *script = PyImport_Import(pName1);
    Py_DECREF(pName1);
    


    //Import all necessary functions: 
    PyObject *getData_func = PyObject_GetAttrString(script, "getData");
    PyObject *printdata_func = PyObject_GetAttrString(script, "printdata");
    PyObject *WriteToDataframe_func = PyObject_GetAttrString(script, "WriteToDataframe");
    PyObject *FOCSinstance_class = PyObject_GetAttrString(FOCS, "FOCSinstance");
    PyObject *FlowNet_class = PyObject_GetAttrString(FOCS, "FlowNet");
    PyObject *FlowOperations_class = PyObject_GetAttrString(FOCS, "FlowOperations");
    PyObject *FOCS_class = PyObject_GetAttrString(FOCS, "FOCS");
   
    
    clock_t t0 = clock();
    //Get data
    struct DataHelp dataToSolve = getDataC(data);
    clock_t t1 = clock();

    //test part///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    InstanceData instance1 = opendata_toC("Data/DEMSdata_FOCS_v1.csv");

    int counter = extract_unique_sorted_times(dataToSolve);

    Graph g(1);
    g.remove_empty();
    
    g.init_focs(instance1, timeStep, instanceSize, randomSample, counter);
    //g.print_graph();   
    g.solve_focs();
    
    clock_t qx = clock();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Convert to Python object to save
    PyObject *dataArgs = PyTuple_Pack(2, dataToSolve.Cols, dataToSolve.Names);
    PyObject *instanceData = PyObject_CallObject(WriteToDataframe_func, dataArgs);
    clock_t t2 = clock();
    //Printdata (dont do this when timing.)
    //PyObject_CallObject(printdata_func, PyTuple_Pack(1, instanceData)); //Print the data obtained
    
    //Build FOCS_instance
    PyObject *focsInstArgs = Py_BuildValue("(Oi)", instanceData, timeStep);
    PyObject *instance = PyObject_CallObject(FOCSinstance_class, focsInstArgs);
    Py_DECREF(focsInstArgs);
    clock_t t3 = clock();
    //Build Flownet
    PyObject *flowNet = PyObject_CallObject(FlowNet_class, NULL);
    clock_t t4 = clock();
    //call FlowNet.focs_instance_to_network
    PyObject *ret = PyObject_CallMethod(flowNet, "focs_instance_to_network", "O", instance);
    clock_t t5 = clock();
    //Build flowOp
    PyObject *graphG = PyObject_GetAttrString(flowNet, "G");
    PyObject *flowOpArgs = Py_BuildValue("(OO)", graphG, instance);
    PyObject *flowOp = PyObject_CallObject(FlowOperations_class, flowOpArgs);
    clock_t t6 = clock();
    empty_flow_func = PyObject_GetAttrString(flowOp, "empty_flow");
    
    //Run FOCS
    PyObject *focsArgs = Py_BuildValue("(OOO)", instance, flowNet, flowOp);
    PyObject *focs = PyObject_CallObject(FOCS_class, focsArgs);
    clock_t t7 = clock();
    //Set MaxFlowAlgorithm to shortest_augmenting_path
    PyObject *nx = PyImport_ImportModule("networkx.algorithms.flow");
    PyObject *sap = PyObject_GetAttrString(nx, "shortest_augmenting_path");
    if (PyObject_SetAttrString(focs, "flow_func", sap) < 0){
        Fail_EXIT("could not set flow_func");
    }

    PyObject *solve_ret = PyObject_CallMethod(focs, "solve_focs", NULL);
    
    clock_t t8 = clock();

    // Give result
    PyObject *objValPy = PyObject_CallMethod(focs, "objective", NULL);
    if (!objValPy) Fail_EXIT("objective() failed");

    double obj_val = PyFloat_AsDouble(objValPy);
    Py_DECREF(objValPy);

    clock_t t9 = clock();

    printf("Objective value = %.9f\n", obj_val);

    //FINISHED

    printf("1) Data loading took            %.5f seconds\n", (double)(t1 - t0) / CLOCKS_PER_SEC);
    printf("X) ALL C PARTS TOOK            %.5f seconds\n", (double)(qx - t1) / CLOCKS_PER_SEC);
    printf("2) Instance creation took       %.5f seconds\n", (double)(t2 - qx) / CLOCKS_PER_SEC);
    printf("3) FlowNet instantiation took   %.5f seconds\n", (double)(t3 - t2) / CLOCKS_PER_SEC);
    printf("4) Building flow network took   %.5f seconds\n", (double)(t4 - t3) / CLOCKS_PER_SEC);
    printf("5) FlowOperations instantiation took %.5f seconds\n", (double)(t5 - t4) / CLOCKS_PER_SEC);
    printf("6) FOCS instantiation took      %.5f seconds\n", (double)(t6 - t5) / CLOCKS_PER_SEC);
    printf("7) Setting flow_func took       %.5f seconds\n", (double)(t7 - t6) / CLOCKS_PER_SEC);
    printf("8) FOCS solve took              %.5f seconds\n", (double)(t8 - t7) / CLOCKS_PER_SEC);
    printf("9) Objective calculation took   %.5f seconds\n", (double)(t9 - t8) / CLOCKS_PER_SEC);

    FILE *out = fopen("timings.csv", "w");
    if (!out) {
        perror("fopen timings.csv");
    } else {
        // header row
        fprintf(out,
            "step,description,seconds\n"
            "1,Data loading,%.5f\n"
            "2,Instance creation,%.5f\n"
            "3,FlowNet instantiation,%.5f\n"
            "4,Building flow network,%.5f\n"
            "5,FlowOperations instantiation,%.5f\n"
            "6,FOCS instantiation,%.5f\n"
            "7,Setting flow_func,%.5f\n"
            "8,FOCS solve,%.5f\n"
            "9,Objective calculation,%.5f\n",
            (double)(t1 - t0) / CLOCKS_PER_SEC,
            (double)(t2 - t1) / CLOCKS_PER_SEC,
            (double)(t3 - t2) / CLOCKS_PER_SEC,
            (double)(t4 - t3) / CLOCKS_PER_SEC,
            (double)(t5 - t4) / CLOCKS_PER_SEC,
            (double)(t6 - t5) / CLOCKS_PER_SEC,
            (double)(t7 - t6) / CLOCKS_PER_SEC,
            (double)(t8 - t7) / CLOCKS_PER_SEC,
            (double)(t9 - t8) / CLOCKS_PER_SEC
        );
        fclose(out);
    }

    //free all things
    Py_XDECREF(flowNet);
    Py_XDECREF(instance);
    Py_XDECREF(ret);
    Py_XDECREF(flowOpArgs);
    Py_XDECREF(graphG);
    Py_XDECREF(flowOp);
    Py_XDECREF(focsArgs);
    Py_XDECREF(sap);
    Py_XDECREF(nx);
    Py_XDECREF(focs);
    Py_XDECREF(dataArgs);

    Py_DECREF(dataToSolve.Cols);
    Py_DECREF(dataToSolve.Names);
    
    Py_XDECREF(empty_flow_func);
    Py_XDECREF(WriteToDataframe_func);
    Py_XDECREF(instanceData);
    Py_XDECREF(FOCS_class);
    Py_XDECREF(FlowOperations_class);
    Py_XDECREF(FlowNet_class);
    Py_XDECREF(FOCSinstance_class);
    Py_XDECREF(printdata_func);
    Py_XDECREF(getData_func);
    Py_XDECREF(script);
    Py_XDECREF(FOCS);
    Py_Finalize();
    //close the csv
    fclose(data);
    return 0;
}

struct DataHelp getDataC(FILE *data) {
        //Get Data in C (this took the longest out of all processes in Python)
        //Ask for filesize
        fseek(data, 0L, SEEK_END);
        long sz = ftell(data);
        // Go back to start
        rewind(data);
        //allocate space
        int line_max = ceil(sz/instanceSize);
        char buffer[line_max];

        PyObject *pyRows  = PyList_New(0);        /* row‑oriented container  */
        PyObject *pyNames = NULL;                 /* column names            */

        int row = 0;
        int ncols = 0;

        while(fgets(buffer, sz, data) && row <= instanceSize) {
            buffer[strcspn(buffer, "\r\n")] = '\0';   /* strip newline */
            
            if (buffer[0] == ',') {
                // shift the whole string one position to the right
                memmove(buffer + strlen("Unnamed: 0"), buffer, line_max - strlen("Unnamed: 0") + 1);  // +1 for the '\0' character

                // Now set the first part of the buffer to "Unnamed: 0,"
                strncpy(buffer, "Unnamed: 0,", strlen("Unnamed: 0,"));
            }      

            char *value = strtok(buffer, ",");     

            if (row == 0) {                       /* ---------- header row -------- */
                pyNames = PyList_New(0);

                while (value) {
                    PyObject *s = PyUnicode_FromString(value);
                    PyList_Append(pyNames, s);
                    Py_DECREF(s);
                    ++ncols;
                    value = strtok(NULL, ",");
                }
            }
            else {                              /* ---------- data row ---------- */
                PyObject *pyRow = PyList_New(0);
                int col = 0;
        
                while (value && col < ncols) {      /* use ncols from header */
                    PyObject *val;
                    char *endptr;

                    // try integer conversion
                    long ival = strtol(value, &endptr, 10);
                    if (*endptr == '\0') {
                        // whole string was an integer
                        val = PyLong_FromLong(ival);
                    } else {
                        // not a pure int, try float
                        double dval = strtod(value, &endptr);
                        if (*endptr == '\0') {
                            val = PyFloat_FromDouble(dval);
                        } else {
                            // fallback to string
                            val = PyUnicode_FromString(value);
                        }
                    }

                    PyList_Append(pyRow, val);
                    Py_DECREF(val);

                    value = strtok(NULL, ",");
                    ++col;
                }
        
                if (col != ncols)
                    Fail_EXIT("row has wrong column count");
        
                PyList_Append(pyRows, pyRow);
                Py_DECREF(pyRow);
            }
            row++;
        }


        Py_ssize_t k = PyList_Size(pyRows);       // number of data rows   
        Py_ssize_t n = PyList_Size(pyNames);      // == ncols               

        PyObject *pyCols = PyList_New(n);         // list of n empty lists  
        for (Py_ssize_t c = 0; c < n; ++c) {
            PyObject *empty = PyList_New(0);
            PyList_SetItem(pyCols, c, empty);     // steals ref
        }

        for (Py_ssize_t r = 0; r < k; ++r) {
            PyObject *rowList = PyList_GetItem(pyRows, r);
            for (Py_ssize_t c = 0; c < n; ++c) {
                PyObject *cell = PyList_GetItem(rowList, c);
                PyObject *colList = PyList_GetItem(pyCols, c);
                PyList_Append(colList, cell);              
            }
        }
        

        
        struct DataHelp a = {pyCols, pyNames};
        Py_DECREF(pyRows);
        return a;
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

int compare_ints(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}

int extract_unique_sorted_times(struct DataHelp data) {

    // Find indices of "t0_timeStep" and "t1_timeStep" columns
    char t0_name[10];
    char t1_name[10];
    snprintf(t0_name, sizeof(t0_name), "t0_%d", timeStep);
    snprintf(t1_name, sizeof(t1_name), "t1_%d", timeStep);

    Py_ssize_t ncols = PyList_Size(data.Names);
    int t0_index = -1, t1_index = -1;
    for (Py_ssize_t i = 0; i < ncols; i++) {    //find the columns which contain t0_timestep and t1_timestep 
        PyObject *name = PyList_GetItem(data.Names, i);  // Borrowed reference
        const char *cname = PyUnicode_AsUTF8(name);
        if (strcmp(cname, t0_name) == 0) {  //if strings equal, we have found the integer.
            t0_index = i;
        } else if (strcmp(cname, t1_name) == 0) {
            t1_index = i;
        }
    }

    if (t0_index == -1 || t1_index == -1) {
        Fail_EXIT("t0_timeStep or t1_timeStep column not found");
    }

    // Get the column lists
    PyObject *t0_col = PyList_GetItem(data.Cols, t0_index);
    PyObject *t1_col = PyList_GetItem(data.Cols, t1_index);

    // Allocate array: 2 * instanceSize ints
    int times[instanceSize*2];
    int counter = 0;

    for (int i = 0; i < instanceSize; i++) {
        PyObject *item0 = PyList_GetItem(t0_col, i);
        PyObject *item1 = PyList_GetItem(t1_col, i);
        times[2 * i] = (int)PyLong_AsLong(item0);   //put the numbers in the 'times' array
        times[2 * i + 1] = (int)PyLong_AsLong(item1);
    }

    // Sort the array
    qsort(times, 2 * instanceSize, sizeof(int), compare_ints);

    FOCS_breakpoints.push_back(times[0]);

    for (int i = 1; i < instanceSize*2; i++) {
        if (times[i] != times[i - 1]) {
            FOCS_breakpoints.push_back(times[i]);;
        }        
    }
    return counter;
}


