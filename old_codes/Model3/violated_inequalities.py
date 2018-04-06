import sys
import time
import re
from read_file_cplex import read_file_dimacs, read_file_formated

# testes if a set of vertices is convergent
def test_if_cover (vertices, list_adj, f):
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

# finds N(S) set neighbors, S list bool selected vertices
def find_neighbors (S, adjacencylist):
    N = len(S)
    neighbor = [0 for _ in range(N)]
    for v in range(N):
        if not S[v]: continue
        for w in adjacencylist[v]:
            neighbor[w] = not S[w]
    list_neighbors = [v for v in range(N) if neighbor[v]]
    return list_neighbors

# finds | N(v) intersect (V/S)|
def neighbors_outside_s (S, v, adjacencylist):
    solution = 0
    for u in adjacencylist[v]:
        if (not S[u]):
            solution += 1
    return solution

# for a given S, tries to fund violated inequalities
def violated_inequality (C0, S, adjacencylist, f):
    # generate every possible S
    # there is a violated inequality if min (c0-1(S) - |N(v) intersection (V/S))
    min_neighbors_outside_S = f[0] + 1
    min_vertex = -1
    C0_S = 0
    # finds number of vertices in the solution and in S
    for v in range(len(f)):
        if S[v]:
            C0_S += C0[v]
    # finds min(f(v) - |N_G(v)) inter (V/S)|
    for v in range(len(f)):
        new_min_val = f[v] - neighbors_outside_s(S, v, adjacencylist)
        if (new_min_val < min_neighbors_outside_S):
            min_neighbors_outside_S = new_min_val
            min_vertex = v
    # a vertex should be choosen as min
    assert(min_vertex != -1)

    if (C0_S < min_neighbors_outside_S):
        # case where we have a violated inequality
        print("Violated: ", min_neighbors_outside_S)
        return min_neighbors_outside_S
    return -1

# Recursively finds all possible vertex combinations of V, and for each of
# this S's tests if it violates any inequality
def generate_S (i, N, S, C0, adjacencylist, f, new_inequalities):
    if i >= N:
        inequality = violated_inequality(C0, S, adjacencylist, f)
        if (inequality != -1):
            new_S = [_ for _ in S]
            new_inequalities.append([new_S, inequality]) #"""
            #print("Inequalities found: ", new_inequalities)
        return
    generate_S(i + 1, N, S, C0, adjacencylist, f, new_inequalities)
    S[i] = True
    generate_S(i + 1, N, S, C0, adjacencylist, f, new_inequalities)
    S[i] = False

# Finds new inequalities for the model, that were violated by the current
# relaxation result
def find_violated_inequalities (C0, adjacencylist, f):
    new_inequalities = []
    S = [False] * len(f)
    generate_S(0, len(f), S, C0, adjacencylist, f, new_inequalities)
    return [len(new_inequalities) != 0, new_inequalities]

if __name__ == "__main__":

    f = [2, 3, 2, 2, 2, 2]
    adjacencylist = [[2, 3, 4, 5], [2, 3, 4, 5], [0, 1], [0, 1], [0, 1], [0, 1]]
    S = [ False ] * len(f)
    #testing find_neighbors
    assert(find_neighbors(S, adjacencylist) == [])
    S[1] = True
    assert(find_neighbors(S, adjacencylist) == [2, 3, 4, 5])
    K = [False] * len(f)
    I = []
    C0 = [False] * len(f)
    C0[0] = True
    C0[1] = True
    C0[2] = True
    generate_S(0, len(f), K, C0, adjacencylist, f, I)
    print(I)
