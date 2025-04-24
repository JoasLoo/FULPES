#include <Python.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int instanceSize = 10; //number of EVs/jobs in instance
int timeStep = 900; //quarterly granularity
//maxFlowAlg = shortest_augmenting_path #alternatively use e.g., edmonds_karp, preflow_push, or dinitz
bool randomSample = true;

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


int main() {
    
    FILE *data; 
    data = fopen("Data/DEMSdata_FOCS_v1.csv","r"); //open data file in C
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
    ////CHECK ALL IMPORTED PYTHON FILES/////
    if (FOCS == NULL || script == NULL) { //If FOCS python has been loaded:
        Fail_EXIT("Couldn't open FOCS.py or script.py");
    } 
    //Get getData function
    PyObject *getData_func = PyObject_GetAttrString(script, "getData");
    //Get printdata function
    PyObject *printdata_func = PyObject_GetAttrString(script, "printdata");
    //Get printdata function
    PyObject *WriteToDataframe_func = PyObject_GetAttrString(script, "WriteToDataframe");
    //Get FOCSinstance class
    PyObject *FOCSinstance_class = PyObject_GetAttrString(FOCS, "FOCSinstance");
    //Get FlowNet class
    PyObject *FlowNet_class = PyObject_GetAttrString(FOCS, "FlowNet");
    //Get FlowOperations
    PyObject *FlowOperations_class = PyObject_GetAttrString(FOCS, "FlowOperations");
    //Get FOCS class
    PyObject *FOCS_class = PyObject_GetAttrString(FOCS, "FOCS");
    //////////////CHECK ALL CLASSES///////////////
    if(FOCS_class == NULL || FlowNet_class == NULL || FOCSinstance_class == NULL || getData_func == NULL) {
        Fail_EXIT("Couldn't read some class or function thing");
    }
    //Get data
    struct DataHelp dataToSolve = getDataC(data);

    //Convert to Python object to save
    PyObject *dataArgs = PyTuple_Pack(2, dataToSolve.Cols, dataToSolve.Names);
    printf("made tuplepack: \n");

    PyObject *instanceData = PyObject_CallObject(WriteToDataframe_func, dataArgs);
    printf("made instancedata: \n");
    //Printdata
    PyObject_CallObject(printdata_func, PyTuple_Pack(1, instanceData)); //Print the data obtained
    
    //Build FOCS_instance
    PyObject *focsInstArgs = Py_BuildValue("(Oi)", instanceData, timeStep);
    PyObject *instance = PyObject_CallObject(FOCSinstance_class, focsInstArgs);
    Py_DECREF(focsInstArgs);
    
    //Build Flownet
    PyObject *flowNet   = PyObject_CallObject(FlowNet_class, NULL);
    
    //call FlowNet.focs_instance_to_network
    PyObject *ret = PyObject_CallMethod(flowNet, "focs_instance_to_network", "O", instance);
    
    //Build flowOp
    PyObject *graphG = PyObject_GetAttrString(flowNet, "G");
    PyObject *flowOpArgs = Py_BuildValue("(OO)", graphG, instance);
    PyObject *flowOp = PyObject_CallObject(FlowOperations_class, flowOpArgs);
    
    //Run FOCS
    PyObject *focsArgs = Py_BuildValue("(OOO)", instance, flowNet, flowOp);
    PyObject *focs = PyObject_CallObject(FOCS_class, focsArgs);
    
    //Set MaxFlowAlgorithm to shortest_augmenting_path
    PyObject *nx = PyImport_ImportModule("networkx.algorithms.flow");
    PyObject *sap = PyObject_GetAttrString(nx, "shortest_augmenting_path");
    if (PyObject_SetAttrString(focs, "flow_func", sap) < 0){
        Fail_EXIT("could not set flow_func");
    }

    
    //Solve FOCS
    //f = focs.solve_focs(MPCstopper=False, MPCcondition=0)

    PyObject *solve_ret = PyObject_CallMethod(focs, "solve_focs", NULL);

    if(!solve_ret) Fail_EXIT("solve_focs failed\n");
    
    // Give result
    PyObject *objValPy = PyObject_CallMethod(focs, "objective", NULL);
    if (!objValPy) Fail_EXIT("objective() failed");

    double obj_val = PyFloat_AsDouble(objValPy);
    Py_DECREF(objValPy);

    printf("Objective value = %.9f\n", obj_val);

    //FINISHED
    
    printf("Run succesful! \n");

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