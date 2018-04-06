"""
preprocessing.py
Before solving an PCI instance using our models, two reduction rules are applied
to the instance, which may generate separated graph components, each component a
new instance.
"""


from read_file_cplex import read_file_dimacs, read_file_formated

REMOVED = -1
INFECTED = -2
DEGREE_ONE = -3
FREE = 0

# This reduction infects vertices that have to be infected in the beginning of
# the process: if N_G(v) < f[v], v is in the solution
def reduction_one(adjacencylist, f, num_neighbors, removed_vertices, type_vertex):
    infected_someone = False
    new_round_need = True
    while (new_round_need):
        new_round_need = False
        for v in range(len(f)):
            if (not removed_vertices[v] and (f[v] <= 0 or
            num_neighbors[v] < f[v])):
                new_round_need = True
                infected_someone = True
                removed_vertices[v] = True
                #if f[v] <= 0, just removed vertex, it is not in the solution
                #print("f ", f[v])
                if f[v] <= 0:
                    type_vertex[v] = REMOVED
                    # else it is in the solution
                else:
                    type_vertex[v] = INFECTED
                for u in adjacencylist[v]:
                    num_neighbors[u] -= 1
                    f[u] -= 1
    return infected_someone

# This reduction removes a vertex if it has 1 neighbor and f[v] = 1
def reduction_two(adjacencylist, f, num_neighbors, removed_vertices, type_vertex):
    removed_vertex = False
    for v in range(len(f)):
        if (not removed_vertices[v] and num_neighbors[v] == 1 and f[v] == 1):
            removed_vertex = True
            removed_vertices[v] = True
            type_vertex[v] = DEGREE_ONE
    return removed_vertex

# Given an instance and a set of vertices to be removed, reduces the
# adjacencylist representation, removing the indicated vertices
def reduced_instance(adjacencylist, f, w, removed_vertices):
    N = len(f)
    adj_matrix = [[False for _ in range(N)] for _ in range(N)]
    for v in range(N):
        for u in adjacencylist[v]:
            adj_matrix[u][v] = True
            adj_matrix[v][u] = True
    new_N = 0
    for v in removed_vertices:
        new_N += 1 - int(v)

    new_vertices_index = [_ for v in range(new_N)]

    new_adjacencylist = [[] for _ in range(new_N)]
    new_f = [_ for _ in range(new_N)]
    new_w = [_ for _ in range(new_N)]
    i = 0
    for v in range(N):
        if (not removed_vertices[v]):
            new_f[i] = f[v]
            new_w[i] = w[v]
            new_vertices_index[i] = v
            assert(f[v] > 0)
            j = 0
            for u in range(N):
                if (not removed_vertices[u]):
                    if adj_matrix[u][v]:
                        new_adjacencylist[i].append(j)
                    j += 1
                    assert(j <= new_N)
            i += 1
            assert(i <= new_N)
    return [new_adjacencylist, new_f, new_w, new_vertices_index]

# Separates the graph in connected components, and its respective instances
def separate_in_connected_instances(adjacencylist, f, w, type_vertex):
    components = []
    visited = [False for v in f]
    for i in range(len(f)):
        if visited[i]: continue
        in_component = [False for v in f]
        queue = [i]
        #print(i)
        while(queue):
            v = queue.pop(0)
            visited[v] = True
            in_component[v] = True
            for u in adjacencylist[v]:
                if not visited[u]:
                    queue.append(u)
                    visited[u] = True
        not_in_component = [(not v) for v in in_component]
        components.append(
            reduced_instance(adjacencylist, f, w, not_in_component))
    return [components, type_vertex]

# detects which instance are looking at, applies the reductions, separate the
# components and returns a list of new instances equivalent to the problem
def data_preprocessing(argv):
    # reads file depending on the format
    if len(argv) == 4:
        f, adjacencylist, w = read_file_dimacs(argv[2], argv[3])
    elif len(argv) == 3:
        f, adjacencylist, w = read_file_formated(argv[2])
    else:
        print("PCI_problem_cplex <model> <datafile> <?option if in dimacs format>")
        exit(0)

    # creates lists that will be used in next functions
    num_neighbors = [(len(adjacencylist[v])) for v in range(len(f))]
    removed_vertices = [False for _ in range(len(f))]
    type_vertex = [FREE for _ in range(len(f))]

    # Use reductions till no vertex is deleted
    vertices_were_removed = True
    while (vertices_were_removed):
        vertices_were_removed = reduction_one(adjacencylist, f, num_neighbors,
            removed_vertices, type_vertex)
        vertices_were_removed = vertices_were_removed or reduction_two(
            adjacencylist, f, num_neighbors, removed_vertices, type_vertex)
    #print("removed ")
    #print(removed_vertices)
    #print(type_vertex)
    for v in range(len(f)):
        assert(removed_vertices[v] == (type_vertex[v] < 0))
    adjacencylist, f, w, new_vertices_index = reduced_instance(adjacencylist, f, w, removed_vertices)
    #print(new_vertices_index)

    return separate_in_connected_instances(adjacencylist, f, w, type_vertex)
