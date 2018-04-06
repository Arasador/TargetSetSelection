#!/usr/bin/python
#
# PCI_problem.py - Model to solve the minimum convergent set
#
# You can run this example at the command line by
# python PCI_problem_cplex.py <filename>

from __future__ import print_function

import sys
import time
import re

import cplex
from cplex.exceptions import CplexSolverError

def get_words(line):
    """Return a list of the tokens in line."""
    line = line.replace("\t", " ")
    line = line.replace("\v", " ")
    line = line.replace("\r", " ")
    line = line.replace("\n", " ")
    while line.count("  "):
        line = line.replace("  ", " ")
    return line.split(" ")

def read_dat_file(filename, option):
    V = 0
    E = 0
    limiarfunction = []
    adjacencylist = []
    adjacencymatrix = []
    continuation = False
    filename_out = "formated_" + option + "_" + filename
    fout = open(filename_out, "w")
    with open(filename) as f:
        for line in f:
            words = get_words(line)
            if (words[0] == "c"):
                continue;
            if (words[0] == "p"):
                fout.write(words[2] + " " + words[3] + "\n")
                V = eval(words[2])
                E = eval(words[3])
                adjacencylist = [[] for _ in range(V)]
                adjacencymatrix = [[False for _ in range(V)] for _ in range(V)]
                counter = 0
                continue;
            words.pop(0);
            u = eval(words[0]) - 1
            v = eval(words[1]) - 1
            adjacencymatrix[u][v] = True
            adjacencymatrix[v][u] = True

    for u in range(V):
        for v in range(V):
            if (adjacencymatrix[u][v]):
                adjacencylist[u].append(v)

    for i in range(V):
        get_limiar_option = {
            "2": 2,
            "degree": max(len(adjacencylist[i]) / 2, 1),
            "e2v": max(E/(2*V), 2)
        }
        k = get_limiar_option.get(option, 2)
        limiarfunction.append(k)
        s = str(k) + "\n"
        fout.write(s)

    for u in range(V):
        for v in range(V):
            if (adjacencymatrix[u][v] and (u < v)):
                fout.write(str(u) + " " + str(v) + "\n")


    fout.close()
    return [limiarfunction, adjacencylist]

def test_if_cover(vertices, list_adj, f):
    V = len(list_adj)
    #print("V ", V, " ", len(vertices))
    neighboors_infected = [0 for _ in range(V)]
    infected = [False for _ in range(V)]
    queue = []
    for v in vertices:
        queue.append(v)
        infected[v] = True
    while queue:
        v = queue.pop(0)
        #print(v)
        for u in list_adj[v]:
            if infected[u]:
                continue
            neighboors_infected[u] += 1
            if neighboors_infected[u] >= f[u]:
                queue.append(u)
                infected[u] = True
    res = True
    for e in infected:
        #print (e, end=", ")
        res = res and e
    #print ()
    return res


def PCI_problem(datafile, option):
    # Read in data file. If no file name is given on the command line
    # we use a default file name. The data we read is
    # limiarfunction   -- a list/array with min value for a vertex be infected
    # adjacencylist    -- a matrix with the neighbors of each vertex
    #

    limiarfunction, adjacencylist = read_dat_file(datafile, option)
    #limiarfunction = [1,2,2,1]
    num_vertices = len(limiarfunction)
    deadline = len(limiarfunction)
    # Create a new (empty) model and populate it below.
    model = cplex.Cplex()

    # Create one binary variable for each vertex/instant pair. The variables
    # model whether vertex is infected in that instant of time or not
    model.variables.add(obj = [1] * num_vertices,
                        lb = [0] * num_vertices,
                        ub = [1] * num_vertices,
                        types = ["B"] * num_vertices)

    for t in range(deadline - 1):
        model.variables.add(obj = [0] * num_vertices,
                            lb = [0] * num_vertices,
                            ub = [1] * num_vertices,
                            types = ["B"] * num_vertices)

    # Create corresponding indices for later use
    x = []
    for t in range(deadline):
        x.append([])
        for v in range(num_vertices):
            x[t].append((t) * (num_vertices) + v)

    # First constraint: The vertex when infected must stay infected
    # x[t-1][v] <= x[t][v] == x[t-1][v] - x[t][v] <= 0
    for t in range(num_vertices - 1):
        for v in range(num_vertices):
            assignment_constraint = cplex.SparsePair(ind=[x[t][v],
                                                     x[t + 1][v]],
                                                     val=[1, -1])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0])


    # Second constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # sum (x[t][u] for u in N(v)) - f(v) <= V * x[t+1][v], for any t
    for t in range(deadline - 1):
        for v in range(num_vertices):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair(ind =v_neighbors + [x[t + 1][v]],
                val=[1] * len(v_neighbors) + [-num_vertices])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[limiarfunction[v]])

    # Third constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # sum (x[t][u] for u in N(v)) - f(v) + 1 <= V * x[t+1][v]
    """for t in range(deadline - 1):
        for v in range(num_vertices):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair(ind = v_neighbors + [x[t + 1][v]],
                 val = [1] * len(v_neighbors) + [-num_vertices])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[limiarfunction[v] - 1])"""

    # Forth constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # + sum (x[t][u] for u in N(v)) + f(v)*x[t][v] >= + f[v] * x[t+1][v]
    # - sum (x[t][u] for u in N(v)) - f(v)*x[t][v] <= - f[v] * x[t+1][v]
    """for t in range(deadline - 1):
        for v in range(num_vertices):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair (ind = v_neighbors + [x[t][v], x[t + 1][v]],
                 val = [-1] * len(v_neighbors) + [-limiarfunction[v], +limiarfunction[v] ])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0]) """


    # Forth constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # + sum (x[t][u] for u in N(v)) + f(v)*x[t][v] >= + f[v] * x[t+1][v]
    # - sum (x[t][u] for u in N(v)) - f(v)*x[t][v] <= - f[v] * x[t+1][v]
    for t in range(deadline - 1):
        for v in range(num_vertices):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair (ind = v_neighbors + [x[t][v], x[t + 1][v]],
                 val = [-1] * len(v_neighbors) + [-num_vertices, +num_vertices ])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[num_vertices - limiarfunction[v]])


    # Fifth constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # + sum (x[t][u] for u in N(v)) = |V|

    assignment_constraint = cplex.SparsePair (ind = [x[deadline - 1][v] for v in range(num_vertices)],
         val = [1] * num_vertices)
    model.linear_constraints.add(lin_expr=[assignment_constraint],
                                 senses=["E"],
                                 rhs=[num_vertices])


    model.parameters.timelimit.set(300.0)
    # Our objective is to minimize cost. Fixed and variable costs
    # have been set when variables were created.
    model.objective.set_sense(model.objective.sense.minimize)
    model.write("mymodel.lp")
    t0 = time.clock()
    model.set_results_stream(None)
    # Solve
    try:
        model.solve()
    except CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        t1 = time.clock()
        solution = model.solution

        # solution.get_status() returns an integer code
        print("Solution status = ", solution.get_status(), ":", end=' ')
        # the following line prints the corresponding string
        print(solution.status[solution.get_status()])

        # Display solution
        out_tables = open("tables/cplex_tables.dat","a")
        out_costs = open("tables/cplex_costs.dat","a")
        out_times = open("tables/cplex_times.dat","a")
        sol_value = solution.get_objective_value()
        print("Total cost = ", sol_value)
        s = datafile + ",\t" + str(sol_value) + ",\t " + str (t1-t0)
        out_tables.write(s)
        s = datafile + ",\t" + str(sol_value)
        out_costs.write(s)
        s = datafile + ",\t" + str(t1 - t0)
        out_costs.write(s)
        print("Time = ", t1 - t0)
        for v in range(num_vertices):
            if (len(adjacencylist[v]) < limiarfunction[v]):
                assert solution.get_values(v) > model.parameters.mip.tolerances.integrality.get()
        #print("Vertices that were initialy infected")
        selected_vertices = []
        for v in range(num_vertices):
            if solution.get_values(v) > model.parameters.mip.tolerances.integrality.get():
                selected_vertices.append(v)
                #print("%d, " % v, end=' ')

        #print()#"""
        print("CPLEX")
        assert(test_if_cover(selected_vertices, adjacencylist, limiarfunction))

if __name__ == "__main__":
    datafile = "../../examples/data/facility.dat"
    if len(sys.argv) != 3:
        print("PCI_problem_cplex <datafile> <option> : ")
        exit(0)
    else:
        datafile = sys.argv[1]
        option = sys.argv[2]
    PCI_problem(datafile, option)
