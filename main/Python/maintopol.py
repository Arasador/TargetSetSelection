from __future__ import print_function
import sys
from TopologicalOrderModel import modeltop
from preprocessing import data_preprocessing
from model1 import PCI_model_1
from model2 import *

INFECTED = -2 # vertex has to be in the final solution
TIMELIMIT = 300

def add_to_table(filename, models_selected, solution, time, num_constraints,
    status, gap, outfile):
    with open(outfile, "a") as myfile:
        myfile.write(str(filename) + ",  ")
        for m in range(len(solution)):
            if models_selected[m] == '1':
                if status[m] == 101:
                    myfile.write("MIP optimal, " )
                elif status[m] == 107:
                    myfile.write("MIP timelimit, " )
                else:
                    myfile.write("MIP ELSE, " )

        for m in range(len(solution)):
            if models_selected[m] == '1':
                myfile.write(str(solution[m]) + ", " )
        for m in range(len(solution)):
            if models_selected[m] == '1':
                myfile.write(str(gap[m]) + ", " )
        for m in range(len(solution)):
            if models_selected[m] == '1':
                myfile.write(str(time[m]) + ", " )
        for m in range(len(solution)):
            if models_selected[m] == '1':
                myfile.write(str(num_constraints[m]) + ", " )
                print("----> Sol[", m, "]   ",  solution[m], " Time ", time[m],
                " Constraints ", num_constraints[m], " GAP ", gap[m])
        myfile.write("\n")


if __name__ == "__main__":
    #apply reduction rules and separate into components
    components, type_vertex = data_preprocessing(sys.argv)
    print("Ended preprocessing")
    infected = 0
    for i in range(len(type_vertex)):
        if type_vertex[i] == INFECTED:
            infected += 1
    #print("Infected ", infected)
    #string of 1's and 0's selects the models we want to run
    models_selected = sys.argv[1]
    models = [# PCI_model_1,
        s_model_initial,
        smallest_s_model,
        smallest_s_weighted,
        dominated_model,
        dominated_model_weighted,
        modeltop ]
    assert(len(models_selected) >= len(models))
    solution, status = [infected] * len(models), [-1] * len(models)
    time, num_constraints = [0] * len(models), [0] * len(models)
    gap = [0] * len(models)

    for component in components:
        adjacencylist = component[0]
        f = component[1]
        w = component[2]
        #print("Component: ", adjacencylist, " f = ", f)
        for m in range(len(models)):
            if models_selected[m] == '1':
                sol, t, num_c, stat, g = models[m](adjacencylist, f, w, TIMELIMIT)
                #print("Solution: ", sol, " t ",  t)
                solution[m] += sol
                time[m] += t
                num_constraints[m] += num_c
                status[m] = max(status[m], stat)
                if status[m] != 101:
                    gap[m] = max(gap[m], g)

    add_to_table(sys.argv[len(sys.argv) - 2], models_selected, solution, time,
        num_constraints, status, gap, "results_PCI")
