from __future__ import print_function

import cplex

import sys

def find_new_adjacency_list(adjacencylist, f, infected, new_N):

    if (new_N == 0):
        return
    N = len(f)

    """ for v in range(N):
        if infected[v]:
            for u in adjacencylist[v]:
                f[u] -= 1
    """
    adj_matrix = [[False for _ in range(N)] for _ in range(N)]
    for v in range(N):
        if infected[v]: continue
        for u in adjacencylist[v]:
            if infected[u]: continue
            adj_matrix[u][v] = True
            adj_matrix[v][u] = True

    new_f = []
    new_adj_list = [[] for _ in range(new_N)]
    new_v = 0
    for v in range(N):
        if infected[v]:
            continue
        assert(f[v] > 0)
        new_f.append(f[v])
        new_u = 0
        for u in range(N):
            if infected[u]: continue
            if adj_matrix[v][u]:
                new_adj_list[new_v].append(new_u)
            new_u += 1
            assert(new_u <= new_N)
        new_v += 1
        assert(new_v <= new_N)
    #adjacencylist = list(new_adj_list)
    #f = list(new_f)
    solution_increase = len(f) - new_N
    return [new_N, new_adj_list, new_f, solution_increase]

# find vertices that are necessarily in the solution
def mandatory_infected_vertices(adjacencylist, f):
    N = len(f)
    deadline = len(f)
    num_neighbors = [len(lis) for lis in adjacencylist]
    initially_infected = [False for _ in range(N)]
    infect = True
    #print(num_neighbors)
    num_infected = 0
    while (infect):
        #print("finding infected")
        infect = False
        for v in range(N):
            if (not initially_infected[v] and (num_neighbors[v] < f[v] or
                f[v] <= 0)):
                num_infected += 1
                initially_infected[v] = True
                infect = True
                for u in adjacencylist[v]:
                    num_neighbors[u] -= 1
                    f[u] -= 1
    #print("num_infected ", num_infected)
    new_N = N - num_infected

    #print(f)
    #print(initially_infected)
    # creates new adjacencylist and f, with equivalent instance without infected
    #vertices
    return find_new_adjacency_list(adjacencylist, f, initially_infected, new_N)

# adds variables that will be used in the program
def add_initial_model_variables(model, N, deadline):
    # Our objective is to minimize cost. Fixed and variable costs
    # have been set when variables were created.
    type_chosen = "S" #"B"
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
    model.variables.add(obj = [0] * N,
                        lb = [0] * N,
                        ub = [1] * N,
                        types = [type_chosen] * N)

# adds problems standard constrainst to the model
def add_problems_basic_constraints (model, x, z, C0, deadline, adjacencylist, f):
    N = len(f)
    # First constraint: The vertex when infected must stay infected
    # x[t-1][v] <= x[t][v] == x[t-1][v] - x[t][v] <= 0
    # was considered unnecessary
    """for t in range(deadline - 1): # so the last t + 1 is for deadline
        for v in range(N):
            assignment_constraint = cplex.SparsePair(ind=[x[t][v], x[t + 1][v]],
                                                     val=[1, -1])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0])
    """

    # Second constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # sum (x[t][u] for u in N(v)) - f(v) + 1 <= N * x[t+1][v], for any t ==
    # sum (x[t][u] for u in N(v)) - N * x[t+1][v] <= + f(v) - 1, for any t
    # may be unnecessary
    """for t in range(deadline - 1):
        for v in range(N):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair(
                ind = v_neighbors + [x[t + 1][v]],
                val = [1] * len(v_neighbors) + [-N])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[f[v] - 1])
    """

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


    # Fifth constraint: at each instant t, the number of infected vertices should
    # be greater or equal then t + C0
    # sum (x[t][u] for u in V) <= t + C0 <- solution previously known
    for t in range(deadline - 1):
            assignment_constraint = cplex.SparsePair (ind = [x[t][v] for v in range(N)],
                 val = [1] * N)
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["G"],
                                         rhs=[t + C0])

    # Sixth constraint: at each instant, one new vertex that was not infected
    # before must be infected in the next instant of time. So, you increase:
    # sum (x[t][u] for u in V) + z[t] <= sum (x[t+1][u] for u in V)
    # or you infected all vertices (notice z is the binary control representing
    # the or operation):
    # sum (x[t][u] for u in V) >= |V|*(1 - z[t])

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
                                     rhs=[N])

# given the violated inequalities found with the relaxation's solution, adds
# them to the model
def add_violated_inequalities (model, inequalities, x, N):
    for inequality in inequalities:
        print("inequality", inequality)
        S = inequality[0]
        value_c0 = inequality[1]
        assignment_constraint = cplex.SparsePair (ind = [x[0][v] for v in range(N) if S[v]],
             val = [1 for v in range(N) if S[v]])
        print(S)
        print("<= ", value_c0)
        model.linear_constraints.add(lin_expr=[assignment_constraint],
                                     senses=["G"],
                                     rhs=[value_c0])
