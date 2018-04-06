""" Models: 1 - simple model first made,
    2 - Using S, constraint directly get after S we got
    3 - Using S, constraint got with smaller S, with random vertices infection
"""



from __future__ import print_function
import sys
from preprocessing import data_preprocessing
from model1 import PCI_model_1
from model2 import PCI_model_2

REMOVED = -1 # vertex automatically infected at the beggining
INFECTED = -2 # vertex has to be in the final solution
DEGREE_ONE = -3 # vertex doesn't make difference in the end result (deg = f = 1)
NUM_MODELS = 3

def add_to_table(filename, solution, time, num_constraints, outfile):
    with open(outfile, "a") as myfile:
        myfile.write(str(filename), end= ",  ")
        for m in range(len(idx)):
            myfile.write(str(solution[m]) + ", " +
            str(time[m]) + ", " + str(num_constraints[m]), end=",  ")



if __name__ == "__main__":
    #apply reduction rules and separate into components
    components, type_vertex = data_preprocessing(sys.argv)

    print("preprocessing infected vertices: ")
    print([i for i in range(len(type_vertex)) if type_vertex[i] == INFECTED])
    print("preprocessing removed vertices: ")
    print([i for i in range(len(type_vertex)) if type_vertex[i] == REMOVED or
        type_vertex[i] == DEGREE_ONE ])

    print(len(components), " components:")
    for c in components:
        print(c)

    num_infected = 0
    for t in type_vertex:
        if t == INFECTED:
            num_infected += 1

    #for each component, solve the problem
    solution = [ num_infected ] * NUM_MODELS
    time = [ 0 ] * NUM_MODELS
    num_constraints = [ 0 ] * NUM_MODELS

    all_model = sys.argv[1] == "all"
    model1 = sys.argv[1] == "m1" or sys.argv[1] == "m12" or sys.argv[1] == "m13" or sys.argv[1] == "m123"
    model2 = sys.argv[1] == "m2" or sys.argv[1] == "m12" or sys.argv[1] == "m23" or sys.argv[1] == "m123"
    model3 = sys.argv[1] == "m3" or sys.argv[1] == "m13" or sys.argv[1] == "m23" or sys.argv[1] == "m123"
    #"""
    for component in components:
        adjacencylist = component[0]
        f = component[1]
        if all_model or model1:
            print("Model 1:")
            sol, t = PCI_model_1(adjacencylist, f)
            solution[0] += sol
            time[0] += t
        if all_model or model2:
            print("Model 2:")
            sol, t, cons = PCI_model_2(adjacencylist, f, "biggerS")
            solution[1] += sol
            time[1] += t
            num_constraints[1] += cons
        if all_model or model3:
            print("Model 3:")
            sol, t, cons = PCI_model_2(adjacencylist, f, "smallerS")
            solution[2] += sol
            time[2] += t
            num_constraints[2] += cons



    if all_model or model1:
        print("----> Sol1 ", solution[0], " Time1 ", time[0])
    if all_model or model2:
        print("----> Sol2 ", solution[1], " Time2 ", time[1], " Constraints ", num_constraints[1])
        #add_to_table( sys.argv[len(argv) - 1], solution, time, num_constraints, "results_PCI")
    if all_model or model3:
        print("----> Sol3 ", solution[2], " Time3 ", time[2], " Constraints ", num_constraints[2])
        #add_to_table(sys.argv[len(argv) - 1], time, num_constraints, "results_PCI")
    #"""
