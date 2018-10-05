#include "PCI_solver.h"
#include "preprocessing.cpp"

#define TIMELIMIT 300


// given solutions, writes results in a file
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

void get_input_from_component(int i, deque<vector<vector<int>>>& components,
  vector<vector<int>>& adjacency_list, vector<int>& f, vector<int>& w) {
  adjacency_list = components[i];
  f = adjacency_list[adjacency_list.size() - 2];
  w = adjacency_list[adjacency_list.size() - 1];
  adjacency_list.pop_back();
  adjacency_list.pop_back();
}

void add_root_relaxation_constraints(PCI_solver& pci_solver) {
  pci_solver.setModelVariables(true);
  //int count = 0;
  while (true) {
    pci_solver.solveProblem();
    if (! pci_solver.recursive_relax()) {
      cout << "break out" << endl;
      break;
    }
    cout << "new constraint added!" << endl;
  }
}


int main (int argc, char** argv) {
  // gets models selected from std input cmd line
  string model_selected = argv[1];
  cout << boolalpha;
  // reads input file and stores it into components vector
  deque<vector<vector<int> > > components = data_preprocessing(argc, argv);
  // possible selected models
  vector<int> objective_values(NUM_MODELS, 0);
  vector<double> times(NUM_MODELS, 0), gaps(NUM_MODELS, 0);

  #ifdef FILE_S_CUTTER_INFO
    ofstream outfile("out_lazyconstraint.txt", ios::app);
    outfile << argv[2] << " \n" ;
    outfile.close();
  #endif
  // list of all models we have
  model models[] = {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
  S_SMALLER_H1, S_SMALLER_H2, S_SMALLER_NEW};

  for (int j = 0; j < components.size(); j ++) {
    vector<vector<int> > adjacency_list;
    vector<int> f , w;
    get_input_from_component(j, components, adjacency_list, f, w);

    for (int i = 0; i < NUM_MODELS; i ++) {
      if (model_selected[i] == '0') continue;
      #ifdef FILE_S_CUTTER_INFO
        ofstream outfile("out_lazyconstraint.txt", ios::app);
        outfile << "\nMODEL " << i << "\n";
        outfile.close();
      #endif
      // builds cplex model
      IloEnv env_pci;
      PCI_solver pci_solver(env_pci, adjacency_list, f, w);
      pci_solver.setModelProblem();
      pci_solver.setCplexSettings(TIMELIMIT);
      //vlDisp, vlEmph, alg, numThreads, vlGap, memory);
      try {
        pci_solver.startAlg(models[i]);
        //one_cut.enforceIntVars();
        vector<vector<int>> lhs;
        vector<int> rhs;
        add_root_relaxation_constraints(pci_solver, lhs, rhs);
        pci_solver.setModelVariables(false);
        pci_solver.setModelProblem ();
        pci_solver.solveProblem();
        pci_solver.endAlg(objective_values[i], times[i], gaps[i]);
        cout << "ended" << endl;
      } catch (IloException& ex) {
        cout << "exception" << endl;
        pci_solver.Texception(objective_values[i], times[i], gaps[i]);
      }
    }
  }
  write_output_file(objective_values, times, gaps, argv[1], argv[2],
    "table_COLORED.dat"); //*/
}
