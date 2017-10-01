# other functions

def test_if_cover(vertices, list_adj, f):
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

def read_file(filename):
    datafile = open(filename, "r")
    line = datafile.readline()
    words = line.split()
    N = eval(words[0])
    M = eval(words[1])
    #reads limiar function from file
    f = []
    for f in range(N):
        line = datafile.readline()
        f.append(eval(line))
        adjacencylist.append([])
    #reads adjacency list from file
    adjacencylist = [[] for _ in range(N)]
    for line in datafile:
        words = line.split()
        u = eval(words[0])
        v = eval(words[1])
        assert(u < v)
        adjacencylist[u].append(v)
        adjacencylist[v].append(u)
    #closes file and return data to model
    datafile.close()
    return [f, adjacencylist]

def find_neighbors(S, adjacencylist):
    N = len(S)
    neighbor = [0 for _ in range(N)]
    for v in range(S):
        if not S[v]: continue
        for w in adjacencylist[v]:
            neighbor[w] = not S[w]
    list_neighbors = [v for v in range(N) if neighbor[v]]
    return list_neighbors

"""def greedy_heuristic(initial_vertex, adjacencylist, f):
    N = len(adjacencylist)
    S = [False for _ in range(N)]
    S[initial_vertex] = True
    resistence_S = f[initial_vertex]
    degree_S = len(adjacencylist[initial_vertex])
    gap_infection = resistence_S - degree_S

    max_w = -1
    max_delta = -1
    neighbors_S = find_neighbors(S, adjacencylist)
    while(True):
        max_w = -1
        # max delta will be maintained
        for w in neighbors_S:
            if (S[w]) continue
            delta_w = 0
            if f[w] < resistence_S:
                delta_w = resistence_S - f[w]
            S[w] = True
            delta_w = delta_w - (len(neighbors_S) - len(find_neighbors(S, adjacencylist)))
            S[w] = False
            if (max_delta < delta_w):
                max_delta = delta_w
                max_w = w
        if max_w <= -1:
            break
        S[max_w] = True
    return [S, max_delta]
#"""

def find_initial_inequalities(adjacencylist, f):
    for initial_v in range(len(adjacencylist)):
        S, delta = greedy_heuristic(initial_v, adjacencylist, f)
        #insert here inequality
