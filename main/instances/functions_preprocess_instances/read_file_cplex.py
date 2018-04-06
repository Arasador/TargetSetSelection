from __future__ import print_function

import sys
import time
import re
import random


def get_words(line):
    """Return a list of the tokens in line."""
    line = line.replace("\t", " ")
    line = line.replace("\v", " ")
    line = line.replace("\r", " ")
    line = line.replace("\n", " ")
    while line.count("  "):
        line = line.replace("  ", " ")
    return line.split(" ")

def read_file_dimacs(filename, option):
    print ("dimacs")
    V = 0
    E = 0
    limiarfunction = []
    adjacencylist = []
    adjacencymatrix = []
    continuation = False
    #filename_out = "formated_" + option + "_" + filename
    #fout = open(filename_out, "w")
    with open(filename) as f:
        for line in f:
            words = get_words(line)
            if (words[0] == "c"):
                continue;
            if (words[0] == "p"):
                #fout.write(words[2] + " " + words[3] + "\n")
                V = eval(words[2])
                E = eval(words[3])
                adjacencylist = [[] for _ in range(V)]
                adjacencymatrix = [[False for _ in range(V)] for _ in range(V)]
                counter = 0
                continue
            words.pop(0)
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
            "degree": max(len(adjacencylist[i]), 1),
            "degree2": max(len(adjacencylist[i]) / 2, 1),
            "e2v": max(E/(2*V), 2)
        }
        k = get_limiar_option.get(option, 2)
        limiarfunction.append(k)
        #s = str(k) + "\n"
        #fout.write(s)

    """for u in range(V):
        for v in range(V):
            if (adjacencymatrix[u][v] and (u < v)):
                fout.write(str(u) + " " + str(v) + "\n")
    """
    #weights = [1] * V
    #fout.close()
    weights = [random.randint(1, 100) for v in range(V)]
    return [limiarfunction, adjacencylist, weights]

def read_file_formated(filename):
    print("formated")
    datafile = open(filename, "r")
    line = datafile.readline()
    words = line.split()
    N = eval(words[0])
    M = eval(words[1])
    #print(N)
    #print(M)
    #reads limiar function from file
    limiarfunction = []
    weights = []
    for i in range(N):
        line = datafile.readline()
        words = line.split()
        limiarfunction.append(eval(words[0]))
        weights.append(eval(words[1]))
        #limiarfunction.append(eval(line))
    #reads adjacency list from file """
    adjacencylist = [[] for _ in range(N)]
    for line in datafile:
        print(line)
        words = line.split()
        u = eval(words[0])
        v = eval(words[1])
        #assert(u < v)
        adjacencylist[u].append(v)
        adjacencylist[v].append(u)
    #closes file and return data to model
    datafile.close()

    return [limiarfunction, adjacencylist, weights]
