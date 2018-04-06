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
import math
from read_file_cplex import read_file_dimacs, read_file_formated
from violated_inequalities import find_violated_inequalities, test_if_cover
from inequalities import *

import cplex
from cplex.exceptions import CplexSolverError


def PCI_problem(N, adjacencylist, f, solution_increase, datafile):
    # Read in data file. If no file name is given on the command line
    # we use a default file name. The data we read is
    # f   -- a list/array with min value for a vertex be infected
    # adjacencylist    -- a matrix with the neighbors of each vertex
    #
    # infects vertices that have to be infected

    print(N)
    print(adjacencylist)
    print(f)


    # variables that control the while loop
    there_are_new_violated_inequalities = True
    violated_inequalities = []
    max_num_iterations = 5
    C0 = 1
    C0_changed = False
    while ((there_are_new_violated_inequalities or C0_changed) and max_num_iterations):
        #print("C0 " ,C0)
        there_are_new_violated_inequalities = False
        max_num_iterations -= 1
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
        #add violated inequalities
        add_violated_inequalities (model, violated_inequalities, x, N)
        #set time limit for convergence
        model.parameters.timelimit.set(300.0)
        #prints the final model for this instance in the file mymodel.lp
        model.write("mymodel.lp")
        #starts counting solver's time
        t0 = time.clock()
        #hides result stream
        model.set_results_stream(None)
        #solve
        try:
            model.solve()
        except CplexSolverError as e:
            print("Exception raised during solve: " + e)
        else:
            #solving termination time
            t1 = time.clock()
            #solution found
            solution = model.solution
            # solution.get_status() returns an integer code
            print("Solution status = ", solution.get_status(), ":", end=' ')
            # the following line prints the corresponding string
            print(solution.status[solution.get_status()])
            # Display solution
            sol_value = solution.get_objective_value()
            out_tables = open("tables/cplex_tables.dat","a")
            s = datafile + ",\t" + str(sol_value) + ",\t " + str (t1-t0)
            out_tables.write(s)
            out_tables.close()
            print("Total cost = ", sol_value, " plus Initially infected" ,
                solution_increase)
            print("Time = ", t1 - t0)
            #for v in range(N):
                #if (len(adjacencylist[v]) < f[v]):
                    #assert solution.get_values(v) > model.parameters.mip.tolerances.integrality.get()
            print("Vertices that were initialy infected")
            C0_var = [False] * N
            for v in range(N):
                print("%f, " %solution.get_values(v), end=' ')
                C0_var[v] = solution.get_values(v)
            #assert(test_if_cover(selected_vertices, adjacencylist, f))
            #find violated inequalities
            there_are_new_violated_inequalities, new_violated_inequalities = \
                find_violated_inequalities(C0_var, adjacencylist, f)
            violated_inequalities += new_violated_inequalities
            #there_are_new_violated_inequalities = False
            #print(violated_inequalities)
            if (not there_are_new_violated_inequalities):
                print("No more violated inequalities")
            #for i in violated_inequalities:
            #    print(i, end=", ")
            #print()
            new_C0 = 0
            for i in range(N):
                new_C0 += C0_var[i]
            new_C0 = int(math.ceil(new_C0))
            C0_changed = new_C0 != C0
            C0 = new_C0
            if (C0_changed):
                print("C0 changed")

def brute_force(i, vertices, adjacencylist, f, min_value):
    if (len(vertices) + 1 >= min_value[0] or i >= len(adjacencylist)) :
        return
    vertices.append(i)
    if (test_if_cover(vertices, adjacencylist, f)) :
        min_value[0] = min(len(vertices), min_value[0])
        vertices.pop()
        return
    vertices.pop()
    brute_force(i + 1, vertices, adjacencylist, f, min_value)
    vertices.append(i)
    brute_force(i + 1, vertices, adjacencylist, f, min_value)
    vertices.pop()

if __name__ == "__main__":
    if len(sys.argv) == 3:
        f, adjacencylist = read_file_dimacs(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        f, adjacencylist = read_file_formated(sys.argv[1])
    else:
        print("PCI_problem_cplex <datafile> <?option if dimacs format>")
        exit(0)
    N, new_adj_list, new_f, solution_increase = \
        mandatory_infected_vertices(adjacencylist, f)
    PCI_problem(N, new_adj_list, new_f, solution_increase, sys.argv[1])
    #v = []
    #min_val = [len(f)]
    #brute_force(0, v, adjacencylist, f, min_val)
    #print("Brute force: ", min_val[0])
