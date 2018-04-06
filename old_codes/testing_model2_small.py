from __future__ import print_function
import sys
from preprocessing import data_preprocessing
from model1 import PCI_model_1
from model2 import PCI_model_2

REMOVED = -1
INFECTED = -2
DEGREE_ONE = -3

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

    #for each component, solve the problem
    solution1, solution2 = 0, 0
    time1, time2 = 0, 0
    model1 = sys.argv[1] == "m1" or sys.argv[1] == "m12"
    model2 = sys.argv[1] == "m2" or sys.argv[1] == "m12"
    #"""
    for component in components:
        adjacencylist = component[0]
        f = component[1]
        if (model1):
            sol_value1, t1 = PCI_model_1(adjacencylist, f)
            solution1 += sol_value1
            time1 += t1
        if (model2):
            sol_value2, t2 = PCI_model_2(adjacencylist, f)
            solution2 += sol_value2
            time2 += t2



    if (model1):
        print("----> Sol1 ", solution1, " Time1 ", time1)
    if (model2):
        print("----> Sol2 ", solution2, " Time2 ", time2)
    #"""
    """
    #s_cutter = S_cutter(adjacencylist, f)
    #infected_vertices = [False, False, False, True, True, False, False, False, False]
    #print(s_cutter.found_new_S_constraints(infected_vertices))

    sol_value, time = PCI_model_2(adjacencylist, f,solution_increase, sys.argv[1])
    out_tables = open("tables/cplex_tables_2.dat","a")
    s = sys.argv[1] + ",\t" + str(sol_value) + ",\t " + str(time) + "\n"
    out_tables.write(s)
    """
