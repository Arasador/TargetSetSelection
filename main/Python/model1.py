from __future__ import print_function

import sys
import time
import re
import math
from preprocessing import data_preprocessing

import cplex
from cplex.exceptions import CplexSolverError,CplexError

def add_initial_model_variables(model, N, deadline):
    # Our objective is to minimize cost. Fixed and variable costs
    # have been set when variables were created.
    type_chosen = "B"
    model.objective.set_sense(model.objective.sense.minimize)
    # variables in first instant of time, used in the objective function
    model.variables.add(obj = [1] * N,
                        lb = [0] * N,
                        ub = [1] * N,
                        types = [type_chosen] * N)


    # variables that do not influence the objective function
    for t in range(1, deadline):
        model.variables.add(obj = [0] * N,
                            lb = [0] * N,
                            ub = [1] * N,
                            types = [type_chosen] * N)

    # variables z, used in inequality seven
    """model.variables.add(obj = [0] * N,
                        lb = [0] * N,
                        ub = [1] * N,
                        types = [type_chosen] * N) """

# adds problems standard constrainst to the model
def add_problems_basic_constraints (model, x, z, C0, deadline, adjacencylist, f):
    N = len(f)
    # First constraint: The vertex when infected must stay infected
    # x[t-1][v] <= x[t][v] == x[t-1][v] - x[t][v] <= 0
    # was considered unnecessary
    for t in range(deadline - 1): # so the last t + 1 is for deadline
        for v in range(N):
            assignment_constraint = cplex.SparsePair(ind=[x[t][v], x[t + 1][v]],
                                                     val=[1, -1])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0])


    # Second constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # sum (x[t][u] for u in N(v)) - f(v) + 1 <= N * x[t+1][v], for any t ==
    # sum (x[t][u] for u in N(v)) - N * x[t+1][v] <= + f(v) - 1, for any t
    # may be unnecessary
    for t in range(deadline - 1):
        for v in range(N):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair(
                ind = v_neighbors + [x[t + 1][v]],
                val = [1] * len(v_neighbors) + [-N])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[f[v] - 1])


    # Third constraint: When a vertex is infected, then it must have the
    #necessary number of neighbors or be infected previously
    # + sum (x[t][u] for u in N(v)) + f(v)*x[t][v] >= + f[v] * x[t+1][v]
    # - sum (x[t][u] for u in N(v)) - f(v)*x[t][v] + f[v] * x[t+1][v] <= 0
    for t in range(deadline - 1):
        for v in range(N):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair (
                ind = v_neighbors + [x[t][v], x[t + 1][v]],
                val = [-1] * len(v_neighbors) + [- f[v], + f[v]])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0])


    # Forth constraint: At the end, all vertices must be infected
    # + sum (x[t][u] for u in N_G(v)) = N
    assignment_constraint = cplex.SparsePair (ind = [x[deadline - 1][v] for v in range(N)],
         val = [1] * N)
    model.linear_constraints.add(lin_expr=[assignment_constraint],
                                 senses=["E"],
                                 rhs=[N])


    # Improvement:
    #Fifth constraint: at each instant t, the number of infected vertices should
    # be greater or equal then t + C0
    # sum (x[t][u] for u in V) <= t + C0 <- solution previously known
    for t in range(deadline - 1):
            assignment_constraint = cplex.SparsePair (ind = [x[t][v] for v in range(N)],
                 val = [1] * N)
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["G"],
                                         rhs=[t + C0])


    # Improvement
    # Sixth constraint: at each instant, one new vertex that was not infected
    # before must be infected in the next instant of time. So, you increase:
    # sum (x[t][u] for u in V) + z[t] <= sum (x[t+1][u] for u in V)
    # or you infected all vertices (notice z is the binary control representing
    # the or operation):
    # sum (x[t][u] for u in V) >= |V|*(1 - z[t])
    """
    for t in range(deadline - 1):
        v_it0 = [x[t][v] for v in range(N)]
        v_it1 = [x[t + 1][v] for v in range(N)]
        assignment_constraint = cplex.SparsePair (
            ind = v_it0 + v_it1 + [z[t]],
            val = [1] * N + [-1] * N + [+1])
        model.linear_constraints.add(lin_expr=[assignment_constraint],
                                     senses=["L"],
                                     rhs=[0])

        assignment_constraint = cplex.SparsePair (
            ind = [x[t][v] for v in range(N)] + [z[t]],
            val = [1] * N + [+ N ])
        model.linear_constraints.add(lin_expr=[assignment_constraint],
                                     senses=["G"],
                                     rhs=[N]) """

# Model 1: basic constraints for PCI problem
def PCI_model_1(adjacencylist, f):
    N = len(f)
    C0 = min(f)
    deadline = N - C0 + 1
    # x and z are the index of the corresponding variable in the model
    x = []
    z = []
    for t in range(deadline):
        x.append([])
        for v in range(N):
            x[t].append((t) * (N) + v)
        z.append(deadline * N + t)

    model = cplex.Cplex()
    #set up the variables used in the model
    add_initial_model_variables(model, N, deadline)
    #add original constraints we have found for this problem
    add_problems_basic_constraints(model, x, z, C0, deadline, adjacencylist, f)
    model.parameters.timelimit.set(3000.0)
    initial_time = model.get_time()
    model.set_results_stream(None)
    try:
        model.solve()
    except CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        #"""
        time = float('%.4f'%(model.get_time() - initial_time))
        solution = model.solution
        sol_value = solution.get_objective_value()
        for i in range(deadline):
            infected_vertices = [(solution.get_values(i * N + v)) for v in range(N)]
            #print(infected_vertices)
        print("\nSolution status = ", solution.get_status(), ":", end=' ')
        # the following line prints the corresponding string
        #print(solution.status[solution.get_status()])

        #print("Total cost = ", sol_value)
        #print("time ", time) #"""
        return [sol_value, time]
