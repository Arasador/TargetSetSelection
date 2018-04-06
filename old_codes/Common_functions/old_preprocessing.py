#MODEL1 PREPROCESSING
#
#
#
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
                f[v] < 0)):
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

def data_preprocessing(argv):
    if len(argv) == 3:
        initial_f, initial_adjacencylist = read_file_dimacs(argv[1], argv[2])
    elif len(argv) == 2:
        initial_f, initial_adjacencylist = read_file_formated(argv[1])
    else:
        print("PCI_problem_cplex <datafile> <?option if dimacs format>")
        exit(0)
    N, adjacencylist, f, solution_increase = \
        mandatory_infected_vertices(initial_adjacencylist, initial_f)
    return [adjacencylist, f, solution_increase]





#MODEL 2 PROPROCESSING

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
                f[v] < 0)):
                num_infected += 1
                initially_infected[v] = True
                #printf("init infev)
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

def data_preprocessing(argv):
    if len(argv) == 3:
        initial_f, initial_adjacencylist = read_file_dimacs(argv[1], argv[2])
    elif len(argv) == 2:
        initial_f, initial_adjacencylist = read_file_formated(argv[1])
    else:
        print("PCI_problem_cplex <datafile> <?option if dimacs format>")
        exit(0)
    f = [len(initial_adjacencylist[v]) for v in range(len(initial_f))]
    N, adjacencylist, f, solution_increase = \
        mandatory_infected_vertices(initial_adjacencylist, initial_f)
    return [adjacencylist, f, solution_increase]


####################################################
"""
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
    return [new_adj_list, new_f, solution_increase]

# find vertices that are necessarily in the solution
def reduce_instance(adjacencylist, f):
    N = len(f)
    deadline = len(f)
    num_neighbors = [len(lis) for lis in adjacencylist]
    initially_infected = [False for _ in range(N)]
    infected_someone = True
    num_infected = 0
    # while finds verteces that can be removed and infected
    while (infected_someone):
        #print("finding infected")
        infected_someone = False
        for v in range(N):
            if (not initially_infected[v] and (num_neighbors[v] < f[v] or
                f[v] < 0)):
                num_infected += 1
                initially_infected[v] = True
                #printf("init infev)
                infected_someone = True
                for u in adjacencylist[v]:
                    num_neighbors[u] -= 1
                    f[u] -= 1

    return find_new_adjacency_list(adjacencylist, f, initially_infected, N - num_infected)
"""
