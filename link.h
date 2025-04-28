#include <Python.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

extern PyObject *empty_flow_func;
extern PyObject *calculate_total_demand_r_func;

void init_FOCS( PyObject *instance ,PyObject *flowNet, PyObject *flowOp, int I_a_count);

void update_network_capacities_g(PyObject *G_r);//,;

void solve_focs_C();