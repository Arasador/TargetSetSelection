from __future__ import print_function

import sys
import time
import re
import math
from preprocessing import data_preprocessing
from s_cutter import S_cutter

import cplex
from cplex.exceptions import CplexSolverError,CplexError
from cplex.callbacks import MIPInfoCallback

class GetInfoCallback(MIPInfoCallback):
    def __call__(self):
        if self.has_incumbent():
            self.gap = min(self.get_MIP_relative_gap(), self.gap)
            #timeused = self.get_time() - self.starttime
            #if timeused > self.timelimit and gap < self.acceptablegap:
            #    print("Good enough solution at", timeused, "sec., gap =",
            #          gap, "%, quitting.")
            #    self.aborted = True
            #    self.abort()

def add_initial_model_variables(model, w, N):
    # Our objective is to minimize cost. Fixed and variable costs
    # have been set when variables were created.
    type_chosen = "B"
    model.objective.set_sense(model.objective.sense.minimize)
    # variables represent vertices in first instant of time, used in obj func
    model.variables.add(obj = w,
                        lb = [0] * N,
                        ub = [1] * N,
                        types = [type_chosen] * N)


    # variables that represent edges, not in obj func
    for v in range(N):
        model.variables.add(obj = [0] * N,
                            lb = [0] * N,
                            ub = [1] * N,
                            types = [type_chosen] * N)

# adds problems standard constrainst to the model
def add_problems_basic_constraints (model, x, e, C0, adjacencylist, f):
    N = len(f)

    #First constraint: sum e_uv >= f[v] * (1 - x_v)
    # == sum e_uv + f[v] * x_v >= f[v]
    for v in range(N):
        v_neighbors = [e[u][v] for u in adjacencylist[v]]
        assignment_constraint = cplex.SparsePair(
            ind = v_neighbors + [x[v]],
            val = [1] * len(v_neighbors) + [f[v]])
        model.linear_constraints.add(lin_expr=[assignment_constraint],
                                     senses=["G"],
                                     rhs=[f[v]])


    #Second constraint: e_uv + e_vu = 1 for every u, v in V distinct
    for u in range(N):
        for v in range(u + 1, N):
            assignment_constraint = cplex.SparsePair(ind=[e[u][v], e[v][u]],
                                                     val=[1, 1])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["E"],
                                         rhs=[1])


    #Third constraint: e_uv + e_vw + e_wu <= 2
    for u in range(N):
        for v in range(N):
            for w in range(N):
                if u == v or u == w or v == w:
                    continue
                assignment_constraint = cplex.SparsePair(
                    ind=[e[u][v], e[v][w], e[w][u]],
                    val=[1, 1, 1])
                model.linear_constraints.add(lin_expr=[assignment_constraint],
                                             senses=["L"],
                                             rhs=[2])

    #improvement: sum x >= C0
    assignment_constraint = cplex.SparsePair(ind = x, val = [1] * N)
    model.linear_constraints.add(lin_expr=[assignment_constraint],
                                 senses=["G"],
                                 rhs=[C0])

# Model 1: basic constraints for PCI problem
def modeltop(adjacencylist, f, w, timelimit):
    N = len(f)
    C0 = min(f)
    deadline = N - C0 + 1
    # x and z are the index of the corresponding variable in the model
    x = []
    e = []
    for u in range(N):
        x.append(u)
        e.append([])
        for v in range(N):
            e[u].append((u + 1) * N + v)
    model = cplex.Cplex()
    #set up the variables used in the model
    add_initial_model_variables(model, w, N)
    #add original constraints we have found for this problem
    add_problems_basic_constraints(model, x, e, C0, adjacencylist, f)

    model.parameters.timelimit.set(timelimit)
    initial_time = model.get_time()
    # get info with mip info callback
    get_info_cb = model.register_callback(GetInfoCallback)
    get_info_cb.gap = 1

    #model.set_results_stream(None)
    try:
        model.solve()
    except CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        #gap = model.getObjValue()
        #print("GAP ", gap)
        time = float('%.4f'%(model.get_time() - initial_time))
        solution = model.solution
        sol_value = solution.get_objective_value()
        infected_vertices = [(solution.get_values(v)) for v in range(N)]
        s_cutter = S_cutter(adjacencylist, f)
        new_f = s_cutter.infect_graph(infected_vertices)
        for i in new_f:
            assert(i < 0)
            #print(infected_vertices)
        print("\nSolution status = ", solution.get_status(), ":", end=' ')
        # the following line prints the corresponding string
        #print(solution.status[solution.get_status()])

        #print("Total cost = ", sol_value)
        #print("time ", time) #"""
        return [sol_value, time, model.linear_constraints.get_num(),
            solution.get_status(), get_info_cb.gap]
