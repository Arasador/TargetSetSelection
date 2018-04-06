#include <vector>
#include <queue>
#include <deque>
#include <stack>
#include <string>
#include <iostream>
#include <fstream>
#include <ilcplex/ilocplex.h>
#include "PCI_solver.h"

#define NUM_MODELS 7
#define TIMELIMIT 300

using namespace std;

// reads data already preprocessed
void read_data (string filename, deque<vector<vector<int> > > &components) {
  int N, M, C, u, v;
  ifstream myfile (filename.c_str());
  if (! myfile.is_open()) {
    cout << "Invalid filename" << endl;
    exit(0);
  }

  myfile >> C;
  for (int c = 0; c < C; c ++) {
    myfile >> N >> M;
    vector<vector<int> > adjacency_list(N);
    vector<int> f(N);
    for (int i = 0; i < N; i ++)
      myfile >> f[i];
    for (int i = 0; i < M; i ++) {
      myfile >> u >> v;
      adjacency_list[u].push_back(v);
      adjacency_list[v].push_back(u);
    }
    adjacency_list.push_back(f);
    components.push_back(adjacency_list);
   }
}

void write_output_file (vector<int> objective_values, vector<double> times,
  vector<double> gaps, string selected_models, string filename, string outname) {
  ofstream outfile(outname, ios::app);
  outfile << filename << ", ";
  for (int i = 0; i < NUM_MODELS; i ++) {
    if (selected_models[i] == '0')
      continue;
    outfile << objective_values[i] << ", " << times[i] << ", " << gaps[i] <<", ";
  }
  outfile << "\n";
  outfile.close();
}

int main (int argc, char** argv) {
  string model_selected = argv[1];
  if (argc != 3 || model_selected.length() < 7) {
    cout << "pci_test <selection of 7 models> <inputfile preprocessed>" << endl;
    exit(0);
  }
  // reads input file and stores it into components vector
  deque<vector<vector<int> > > components;
  read_data(argv[2], components);
  // possible selected models
  int models[] = {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
    S_SMALLER_H1, S_SMALLER_H2};
  vector<int> objective_values(NUM_MODELS, 0);
  vector<double> times(NUM_MODELS, 0), gaps(NUM_MODELS, 0);
  for (int i = 0; i < components.size(); i ++) {
    // for each component runs the model
    //print_matrix(components[i], "Components: \n");
    // reads data from component vector
    vector<vector<int> > adjacency_list = components[i];
    vector<int> f  = adjacency_list[adjacency_list.size() - 1];
    adjacency_list.pop_back();
    components.pop_front();
    for (int i = 0; i < NUM_MODELS; i ++) {
      if (model_selected[i] == '0')
        continue;
      IloEnv env_pci;
      PCI_solver pci_solver(env_pci, adjacency_list, f);
      pci_solver.setCplexSettings(TIMELIMIT);//vlDisp, vlEmph, alg, numThreads, vlGap, memory);
      pci_solver.setModelProblem();
      try {
        pci_solver.startAlg(models[i]);
        //one_cut.enforceIntVars();
        pci_solver.solveProblem();
        pci_solver.endAlg(objective_values[i], times[i], gaps[i]);
      } catch (IloException& ex) {
        pci_solver.Texception();
      }
    }
  }
  write_output_file(objective_values, times, gaps, argv[1], argv[2],
    "output_models.dat");
}
