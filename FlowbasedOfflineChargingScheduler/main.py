#    Flow-based Offline Charging Scheduler (FOCS)
#    For scheduling of large groups of EVs
#    Part of SmoothEMS met GridShield project - code developed by 
#    Leoni Winschermann, University of Twente, l.winschermann@utwente.nl
#    
#    Copyright (C) 2024 CAES and MOR Groups, University of Twente, Enschede, The Netherlands
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
#    USA
import time
import networkx as nx
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

import random

from LP import LP
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance
from Bookkeeping import Bookkeeping

repetitions = 2 #Number of iterations to find average time.

total_load = 0
total_graph = 0
total_init = 0
total_solve = 0



instanceSize = 600 #number of EVs/jobs in instance
timeStep = 900 #quarterly granularity
maxFlowAlg = edmonds_karp #alternatively use e.g., edmonds_karp, preflow_push, or dinitz
randomSample = True
for i in range(repetitions):
    t0 = time.perf_counter()
    # Real Training data
    instanceData = pd.read_csv('../Data/ev_session_data_OR.csv') #ev_session_data_OR  DEMSdata_FOCS_v1

    t1 = time.perf_counter()

    if not randomSample:
        instance = FOCSinstance(instanceData[:instanceSize], timeStep)
    if randomSample:
        sample = sorted(random.sample(range(0,len(instanceData)), instanceSize))
        instance = FOCSinstance(instanceData.iloc[sample], timeStep)  
    t2 = time.perf_counter()
    #print('How instance is formatted:\n', instanceData.iloc[sample])
    '''--------------start FOCS--------------'''
    flowNet = FlowNet()
    flowNet.focs_instance_to_network(instance)
    flowOp = FlowOperations(flowNet.G, instance)
    t3 = time.perf_counter()
    focs = FOCS(instance, flowNet, flowOp)
    focs.flow_func = maxFlowAlg
    t4 = time.perf_counter()
    f = focs.solve_focs(MPCstopper=False, MPCcondition=0)
    t5 = time.perf_counter()

    obj_val = focs.objective()
    #print('FOCS objective value = ', obj_val)

    total_load += (t2 - t1) * 1_000_000
    total_graph += (t3 - t2) * 1_000_000
    total_init  += (t4 - t3) * 1_000_000
    total_solve += (t5 - t4) * 1_000_000


#print('FOCS flow (schedule in middle edge layer): \n', focs.f)

print(f"Avg load time (us): {total_load / repetitions}")
print(f"Avg graph setup time (us): {total_graph / repetitions}")
print(f"Avg init time (us): {total_init / repetitions}")
print(f"Avg solve time (us): {total_solve / repetitions}")



#optional operations
#instance.toy_instance_2() #overwrites (potentially empty) FOCS instance with toy instance
#instance.reduced_problem(f, start_h = 12) #creates instance remaining at noon
#instance.validate_solution(f)
#flowOp.validate_global_cap()
#focsMPC.validate_MPC()
#lp.print_results()
#focs.write()


'''
TODO
General
    synthetic instances
Gurobi
FOCS
    warm-start function
'''