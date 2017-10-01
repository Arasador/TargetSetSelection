from __future__ import print_function

import cplex

import sys

def mandatory_infected_vertices(adjacencylist, f):
    N = len(f)
    deadline = len(f)

    num_neighbors = [len(lis) for lis in adjacencylist]
    initially_infected = [False for _ in range(N)]
    infect = True
    while (infect):
        print("finding infected")
        infect = False
        for v in range(N):
            if (not initially_infected[v] and num_neighbors[v] < f[v]):
                initially_infected[v] = True
                infect = True
                for u in adjacencylist[v]:
                    num_neighbors[u] -= 1
    return initially_infected

def add_initial_model_variables(model, N, deadline):
    # Our objective is to minimize cost. Fixed and variable costs
    # have been set when variables were created.
    model.objective.set_sense(model.objective.sense.minimize)
    # variables in first instant of time, used in the objective function
    model.variables.add(obj = [1] * N,
                        lb = [0] * N,
                        ub = [1] * N,
                        types = ["B"] * N)


    # variables that do not influence the objective function
    for t in range(1, deadline):
        model.variables.add(obj = [0] * N,
                            lb = [0] * N,
                            ub = [1] * N,
                            types = ["S"] * N)

    # variables z, used in inequality seven
    model.variables.add(obj = [0] * N,
                        lb = [0] * N,
                        ub = [1] * N,
                        types = ["S"] * N)

def add_problems_basic_constraints (model, x, z, C0, adjacencylist, f, initially_infected):
    N = len(f)
    deadline = len(f) - C0

    # First constraint: The vertex when infected must stay infected
    # x[t-1][v] <= x[t][v] == x[t-1][v] - x[t][v] <= 0
    # was considered unnecessary
    """for t in range(deadline - 1):
        for v in range(N):
            assignment_constraint = cplex.SparsePair(ind=[x[t][v],
                                                     x[t + 1][v]],
                                                     val=[1, -1])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0])
    """

    # Second constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # sum (x[t][u] for u in N(v)) - f(v) <= V * x[t+1][v], for any t
    # may be unnecessary
    for t in range(deadline - 1):
        for v in range(N):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair(ind =v_neighbors + [x[t + 1][v]],
                val=[1] * len(v_neighbors) + [-N])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[f[v]])

    # Third constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # sum (x[t][u] for u in N(v)) - f(v) + 1 <= V * x[t+1][v]
    """for t in range(deadline - 1):
        for v in range(N):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair(ind = v_neighbors + [x[t + 1][v]],
                 val = [1] * len(v_neighbors) + [-N])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[f[v] - 1])"""

    # Forth constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # + sum (x[t][u] for u in N(v)) + f(v)*x[t][v] >= + f[v] * x[t+1][v]
    # - sum (x[t][u] for u in N(v)) - f(v)*x[t][v] <= - f[v] * x[t+1][v]
    """for t in range(deadline - 1):
        for v in range(N):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair (ind = v_neighbors + [x[t][v], x[t + 1][v]],
                 val = [-1] * len(v_neighbors) + [-f[v], +f[v] ])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0]) """


    # Forth constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # + sum (x[t][u] for u in N(v)) + f(v)*x[t][v] >= + f[v] * x[t+1][v]
    # - sum (x[t][u] for u in N(v)) - f(v)*x[t][v] <= - f[v] * x[t+1][v]
    for t in range(deadline - 1):
        for v in range(N):
            v_neighbors = [x[t][u] for u in adjacencylist[v]]
            assignment_constraint = cplex.SparsePair (ind = v_neighbors + [x[t][v], x[t + 1][v]],
                 val = [-1] * len(v_neighbors) + [-f[v], +f[v] ])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["L"],
                                         rhs=[0])


    # Fifth constraint: When a vertex has the necessary amount of infected
    # neighbors, then it must get infected
    # + sum (x[t][u] for u in N(v)) = |V|

    assignment_constraint = cplex.SparsePair (ind = [x[deadline - 1][v] for v in range(N)],
         val = [1] * N)
    model.linear_constraints.add(lin_expr=[assignment_constraint],
                                 senses=["E"],
                                 rhs=[N])


    # Sixth constraint: at each instant, one new vertex must be infected
    #
    # sum (x[t][u] for u in V) <= t + 1
    for t in range(deadline - 1):
            assignment_constraint = cplex.SparsePair (ind = [x[t][v] for v in range(N)],
                 val = [1] * N)
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["G"],
                                         rhs=[t + C0])

    # Seventh constraint: at each instant, one new vertex that was not infected
    # before must be infected
    # sum (x[t][u] for u in V) + z[t] <= sum (x[t+1][u] for u in V)
    # sum (x[t][u] for u in V) >= |V|*(1 - z[t])

    for t in range(deadline - 1):
        v_it0 = [x[t][v] for v in range(N)]
        v_it1 = [x[t + 1][v] for v in range(N)]
        assignment_constraint = cplex.SparsePair (ind = v_it0 + v_it1 + [z[t]],
             val = [1] * N + [-1] * N + [+1])
        model.linear_constraints.add(lin_expr=[assignment_constraint],
                                     senses=["L"],
                                     rhs=[0])

        assignment_constraint = cplex.SparsePair (ind = [x[t][v] for v in range(N)] + [z[t]],
             val = [1] * N + [+ N ])
        model.linear_constraints.add(lin_expr=[assignment_constraint],
                                     senses=["G"],
                                     rhs=[N])

    # initially_infected must be infected
    print("Initially infected: ", end ="")
    count_intially_infected = 0
    for v in range(N):
        if (initially_infected[v]):
            print(v, end=", ")
            count_intially_infected += 1
            assignment_constraint = cplex.SparsePair (ind = [x[0][v]],
                 val = [1])
            model.linear_constraints.add(lin_expr=[assignment_constraint],
                                         senses=["E"],
                                         rhs=[1])
    print("\n Total: ", count_intially_infected)

#def update_inequality_round_increase (model, C0, x, adjacencylist, f):


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
