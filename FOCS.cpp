#include "link.h"
using namespace std;
//Creating some global variables to save a FOCS instance.
PyObject *instance;
PyObject *flowNet; 
PyObject *flowOp;
vector<int> I_p;    //temp array to store where we are. If empty, done.
vector<int> I_a;    //array counting from 0 to N, where N is the length of: sorted non-double values from t0_<timestep> and t1_<timestep> columns
bool global_cap_active; 

vector<int> I_crit;
char flow_func[] = "shortest_augmenting_path";

PyObject *G;    //initialized from python
PyObject *G_r;  //initialized from python
PyObject *G_rk;
PyObject *total_demand_r;
double *flow_val;
PyObject *f;
//CSR
//https://www.youtube.com/watch?v=tq3zPnrQIpU What are graphs
//https://www.youtube.com/watch?v=nvRkFi8rbOM graphs in C and C++ (replacement for G, G_r and G_rk)
//We have about a 4% coverage of edges (meaning that if we were to make an Adjacency matrix, only 4% of it would be fully useful), therefore we use an Adjacency list.
//For each list item we want to have the int FROM, int TO, int capacity, int usage? More??
//Functions: Init, add_flow, init_focs? (??), reduce_network (??), length_sum_intervals (??), partial_flow (??), calculate_total_demand_r (??).
// all with (??) I have to figure out what it does and why.
//For faster reading: Maybe a sort function? so that the list is sorted based on the FROM or TO int?
//TODO: Make a struct for the digraph
//TODO: Make a vector of the struct
//TODO: Create the functions that are known
//TODO: Find out what the unknown functions should do
//TODO: Make unknown functions.
//TODO: Find out what instance is about.
//TODO: Implement instance in C

//When all these TODO's are done, the code should be sped up sufficiently, then try to implement multithreading.

bool terminate = false;
int rd; //Round counter
int it; //iteration counter

//Variables for solve_focs
double err = 0.0000001; 
bool MPCstopper = false;
int MPCcondition = 0;


void update_I_crit();

    void init_FOCS( PyObject *instance_obtained, PyObject *flowNet_obtained, PyObject *flowOp_obtained, int I_a_count) {
        rd = 0; 
        it = 0;
        terminate = false;
        instance = instance_obtained;
        flowNet = flowNet_obtained;
        flowOp = flowOp_obtained;
        //no need to init space for I_a, I_p and I_crit because theyre global vectors.
        for (int j = 0; j < I_a_count - 1; j++) {
            I_a[j] = j;
        }
        global_cap_active = false;
        G = PyObject_GetAttrString(flowNet, "G");
        G_r = G;
        //total_demand_r = PyObject_CallObject(calculate_total_demand_r_func, PyTuple_Pack(1, G_r));
        flow_val = 0;
        f = PyObject_CallObject(empty_flow_func, PyTuple_Pack(1, G));
    }



    /*#update network capacities (g function in paper)
    def update_network_capacities_g(self,G_r,flow_val): #doesn't have to be same G_r as in self. May be G_rk
        if self.it > 0:
            #print("update capacities of round ", self.rd, "iteration ", self.it)
            demand = self.maxDiff
            demand_normalized = demand/self.flowOp.length_sum_intervals(self.I_a,self.instance.len_i)
            for i in self.I_a:
                G_r["i{}".format(i)]["t"]["capacity"] += demand_normalized * self.instance.len_i[i]
        else:
            #print("initialize capacities of round ", self.rd)
            demand = sum([G_r["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])
            demand_normalized = demand/self.flowOp.length_sum_intervals(self.I_a, self.instance.len_i)         
            for i in self.I_a:
                G_r["i{}".format(i)]["t"]["capacity"] = demand_normalized * self.instance.len_i[i]

        #if global maximum power constraint is active, trigger it here
        if self.instance.global_cap_active:
            for i in self.I_a:
                G_r["i{}".format(i)]["t"]["capacity"] = min(G_r["i{}".format(i)]["t"]["capacity"], self.instance.global_cap[i]*self.instance.len_i[i]*self.instance.tau/self.instance.timeStep)        
        return G_r
    */

    void update_network_capacities_g(PyObject *G_r) {
        
    }

    /*
    def solve_focs(self, err = 0.0000001, MPCstopper = False, MPCcondition = 0):
        while not self.terminate:
            #begin round
            #initiate capacities
            #print("Still working \n")
            if self.it == 0:
                G_rk = copy.deepcopy(self.G_r)
            G_rk = self.update_network_capacities_g(G_rk, self.flow_val)

            #determine max flow
            self.flow_val, flow_dict = nx.maximum_flow(G_rk, "s", "t", flow_func=self.flow_func)
            self.maxDiff = sum([G_rk["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs]) - self.flow_val
            if self.total_demand_r - self.flow_val < err:
                #end round

                #Check for subcrit. For some instances, there are still active intervals that only now don't reach the max anymore #FIXME (for now)
                subCrit_mask = [G_rk["i{}".format(i)]["t"]["capacity"]-flow_dict["i{}".format(i)]["t"] > err for i in self.I_a]
                subCrit = [self.I_a[i] for i in range(0,len(self.I_a)) if subCrit_mask[i]]
                
                #Update I_p
                self.I_p += subCrit
                #Update I_crit
                self.I_crit = self.I_crit + [sorted([self.I_a[i] for i in range(0,len(self.I_a)) if not subCrit_mask[i]])]

                #check terminate
                if len(self.I_p) == 0:
                    self.terminate = True
                    #Add critical flow to current flow
                    self.f = self.flowOp.add_flows(self.f, flow_dict, G_rk)
                else:
                    #determine critical flow
                    I_crit_r = self.I_crit[-1] 
                    flow_dict_crit = self.flowOp.partial_flow(flow_dict, I_crit_r, self.I_a, self.instance.jobs, self.instance.J) #FIXME check if J needs to be reduced?
                    
                    #Add critical flow to current flow
                    self.f = self.flowOp.add_flows(self.f,flow_dict_crit,G_rk)

                    #Reduce flow network by critical flow
                    self.G_r = self.flowOp.reduce_network(self.G_r,flow_dict_crit,I_crit_r)

                    #Update active and parked sets
                    self.I_a = sorted(self.I_p)
                    self.I_p = []

                    if MPCstopper:
                        self.MPCcondition = MPCcondition
                        if MPCcondition < min(self.I_a):
                            # print('[MESSAGE]: Solution optimal for first {} timesteps only.'.format(MPCcondition + 1))
                            return self.f

                    #amount of work yet to be scheduled
                    self.total_demand_r = sum([self.G_r["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])

                    self.rd += 1
                    self.it = 0
            else:
                #initiate next iteration
                #Determine subcritical sets
                if self.instance.global_cap_active:
                    subCrit_mask = [(G_rk["i{}".format(i)]["t"]["capacity"] - flow_dict["i{}".format(i)]["t"] > err) or (self.instance.global_cap[i]*self.instance.len_i[i]*self.instance.tau/self.instance.timeStep == G_rk["i{}".format(i)]["t"]["capacity"]) for i in self.I_a]
                else:
                    subCrit_mask = [G_rk["i{}".format(i)]["t"]["capacity"] - flow_dict["i{}".format(i)]["t"] > err for i in self.I_a]
                subCrit = [self.I_a[i] for i in range(0,len(self.I_a)) if subCrit_mask[i]]

                #Reduce network G_rk
                flow_dict_sub = self.flowOp.partial_flow(flow_dict, subCrit, self.I_a, self.instance.jobs, self.instance.J)
                G_rk = self.flowOp.reduce_network(G_rk,flow_dict_sub,subCrit)
                self.total_demand_r = sum([G_rk["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])
                
                #Update I_p
                self.I_p += subCrit

                #Update I_a
                self.I_a = sorted([self.I_a[i] for i in range(0,len(self.I_a)) if not subCrit_mask[i]])
                self.it += 1

        return self.f*/


        void solve_focs_C() {
            while (!terminate) {
                // begin round
                if(it == 0) {
                    G_rk = G_r;
                }
                update_network_capacities_g(G_rk);

                // determine max flow

                if(PyLong_AsDouble(total_demand_r)-*flow_val < err) {
                    // end the round
                    if (I_p.size() == 0) {
                        // terminate
                        terminate = true;
                        //f = flowOp.addflows? //in python first, later will make in C++
                    }
                    else {
                        // determine critical flow
                    }
                }
                else {
                    // initiate next iteration
                }
                terminate = true;
            }
        }