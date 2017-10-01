# PCI_problem_cplex <data_file> <limiar_option>
# Limiar options
# "2": 2,
# "degree": max(len(adjacencylist[i]) / 2, 1),
# "e2v": max(E/(2*V), 2)

def get_words(line):
    """Return a list of the tokens in line."""
    line = line.replace("\t", " ")
    line = line.replace("\v", " ")
    line = line.replace("\r", " ")
    line = line.replace("\n", " ")
    while line.count("  "):
        line = line.replace("  ", " ")
    return line.split(" ")

def read_dat_file(filename, option):
    V = 0   E = 0
    limiarfunction = []
    adjacencylist = []
    adjacencymatrix = []
    continuation = False
    filename_out = "formated_" + option + "_" + filename
    fout = open(filename_out, "w")
    with open(filename) as f:
        for line in f:
            words = get_words(line)
            if (words[0] == "c"):
                continue;
            if (words[0] == "p"):
                fout.write(words[2] + " " + words[3] + "\n")
                V = eval(words[2])
                E = eval(words[3])
                adjacencylist = [[] for _ in range(V)]
                adjacencymatrix = [[False for _ in range(V)] for _ in range(V)]
                counter = 0
                continue;
            words.pop(0);
            u = eval(words[0]) - 1
            v = eval(words[1]) - 1
            adjacencymatrix[u][v] = True
            adjacencymatrix[v][u] = True

    for u in range(V):
        for v in range(u + 1, V):
            if (adjacencymatrix[u][v]):
                adjacencylist[u].append(v)
                adjacencylist[v].append(u)

    for i in range(V):
        get_limiar_option = {
            "2": 2,
            "degree": max(len(adjacencylist[i]) / 2, 1),
            "e2v": max(E/(2*V), 2)
        }
        k = get_limiar_option.get(option, 2)
        limiarfunction.append(k)
        s = str(k) + "\n"
        fout.write(s)

    for u in range(V):
        for v in range(u + 1, V):
            if (adjacencymatrix[u][v] and (u < v)):
                fout.write(str(u) + " " + str(v) + "\n")
    
    fout.close()

if __name__ == "__main__":
    datafile = "../../examples/data/facility.dat"
    if len(sys.argv) != 3:
        print("PCI_problem_cplex <data_file> <limiar_option> : ")
        exit(0)
    else:
        datafile = sys.argv[1]
        option = sys.argv[2]
        read_dat_file(datafile, option)
