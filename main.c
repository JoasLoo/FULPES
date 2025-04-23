/*import networkx as nx
import copy
from networkx.algorithms.flow import shortest_augmenting_path
from networkx.algorithms.flow import edmonds_karp
from networkx.algorithms.flow import preflow_push
from networkx.algorithms.flow import dinitz
import math
import matplotlib.pyplot as plt
import pandas as pd
import statistics
import csv

#for gurobi model
from gurobipy import *
from itertools import product

import os, sys
import datetime
import time
import random

from LP import LP
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance
from Bookkeeping import Bookkeeping

instanceSize = 10 #number of EVs/jobs in instance
timeStep = 900 #quarterly granularity
maxFlowAlg = shortest_augmenting_path #alternatively use e.g., edmonds_karp, preflow_push, or dinitz
randomSample = True

# Real Training data
instanceData = pd.read_excel('../Data/ev_session_data_OR.xlsx')

if not randomSample:
    instance = FOCSinstance(instanceData[:instanceSize], timeStep)
if randomSample:
    sample = sorted(random.sample(range(0,len(instanceData)), instanceSize))
    instance = FOCSinstance(instanceData.iloc[sample], timeStep)  

'''--------------start FOCS--------------'''
flowNet = FlowNet()
flowNet.focs_instance_to_network(instance)
flowOp = FlowOperations(flowNet.G, instance)
focs = FOCS(instance, flowNet, flowOp)
focs.flow_func = shortest_augmenting_path
f = focs.solve_focs(MPCstopper=False, MPCcondition=0)

obj_val = focs.objective()

print('FOCS objective value = ', obj_val)
print('FOCS flow (schedule in middle edge layer): \n', focs.f)*/

#include <Python.h>
#include <stdio.h>
#include <stdbool.h>

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

int main() {
    Py_Initialize();
    PyObject *pName = PyUnicode_DecodeFSDefault("FOCS");
    PyObject *FOCS = PyImport_Import(pName);     //Import "FOCS".py
    Py_DECREF(pName);
    ////CHECK/////
    if (FOCS == NULL) { //If FOCS python has been loaded:
        Fail_EXIT("Couldn't open FOCS.py");
    } 

    PyObject *FOCS_class = PyObject_GetAttrString(FOCS, "FOCS"); //get hello func

    ////CHECK/////
    if(FOCS_class == NULL) {
        Fail_EXIT("Couldn't read FOCS_class");
    }

    

    printf("Run succesful! \n");
    Py_DECREF(FOCS_class);
    Py_DECREF(FOCS);
    Py_Finalize();
    return 0;
}





////////////////////////////////////////////////example code
/*int main() {
    Py_Initialize(); //Initialize the Python interpreter.

    // Load script.py
    PyObject *pName = PyUnicode_DecodeFSDefault("script");
    PyObject *pModule = PyImport_Import(pName);     //Import "script".py
    Py_DECREF(pName);                               //Free memory used by pName.

    if (pModule != NULL) {
        // Call say_hello()
        PyObject *pFunc_hello = PyObject_GetAttrString(pModule, "say_hello"); //get hello func
        PyObject *pFunc_add = PyObject_GetAttrString(pModule, "add"); //get add function from script.py
        PyObject *args = PyTuple_Pack(3, PyLong_FromLong(3), PyLong_FromLong(5), PyLong_FromLong(100));

        //if (PyCallable_Check(pFunc_hello)) {                                  //check if the func is gotten properly
            PyObject_CallObject(pFunc_hello, NULL);                           //call function with no input args
        //}
        
        if (PyCallable_Check(pFunc_add)) {
            PyObject *pValue = PyObject_CallObject(pFunc_add, args);
            //
            if (pValue != NULL) {
                long result = PyLong_AsLong(pValue);
                printf("Result of add(3, 5): %ld\n", result);
                PyTuple_SetItem(args, 2, PyLong_FromLong(result)); //start counting from 0, so we change the third variable
                Py_DECREF(pValue);
            }
        }
        
        for (int i = 0; i < 0; i++) {
            long result1 = PyLong_AsLong(PyObject_CallObject(pFunc_add, args));
            printf("Result of add(3, 5, 100): %ld\n", result1);
            PyTuple_SetItem(args, 2, PyLong_FromLong(result1));
        }

        //PyObject_CallObject(pFunc_hello, NULL); 
        //free up space
        Py_XDECREF(pFunc_hello);
        Py_XDECREF(pFunc_add);
        Py_XDECREF(args);
        Py_XDECREF(pModule);
    } else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"script.py\"\n");
    }

    Py_Finalize();  //exit the python environment
    return 0;
}*/
