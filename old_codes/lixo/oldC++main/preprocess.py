from preprocessing import data_preprocessing
import sys

if __name__ == "__main__":
    components, type_vertex = data_preprocessing(sys.argv)
    with open("preprocessed_instances/" + str(sys.argv[2]), "w+") as myfile:
        myfile.write(str(len(components)) + "\n")
        for component in components:
            adjacency_list = component[0]
            f = component[1]
            w = component[2]
            s = ""
            M = 0
            for f_i in f:
                s += str(f_i) + "\n"
            for v in range(len(f)):
                for u in adjacency_list[v]:
                    s += str(v) + " " + str(u) + "\n"
                    M += 1
            myfile.write(str(len(f)) + " " + str(M) + "\n" + s)