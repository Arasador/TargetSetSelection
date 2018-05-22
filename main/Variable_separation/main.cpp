#include "preprocessing.cpp"
#include <iostream>
#include <deque>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  deque<vector<vector<int> > > components = data_preprocessing(argc, argv);
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
  for (int i = 0; i < components.size(); i ++) {

    vector<vector<int> > adjacency_list = components[i];
    vector<int> f = adjacency_list[adjacency_list.size() - 2];
    vector<int> w = adjacency_list[adjacency_list.size() - 1];
    adjacency_list.pop_back();
    adjacency_list.pop_back();

    for (int i = 0; i < NUM_MODELS; i ++) {
      if (model_selected[i] == '0')
        continue;

      #ifdef FILE_S_CUTTER_INFO
        ofstream outfile("out_lazyconstraint.txt", ios::app);
        outfile << "\nMODEL " << i << "\n";
        outfile.close();
      #endif

      IloEnv env_pci;
      PCI_solver pci_solver(env_pci, adjacency_list, f, w);
      pci_solver.setModelProblem();
      pci_solver.setCplexSettings(TIMELIMIT);
      //vlDisp, vlEmph, alg, numThreads, vlGap, memory);
      try {
        pci_solver.startAlg(models[i]);
        //one_cut.enforceIntVars();
        pci_solver.solveProblem();
        pci_solver.endAlg(objective_values[i], times[i], gaps[i]);
      } catch (IloException& ex) {
        pci_solver.Texception(objective_values[i], times[i], gaps[i]);
      }
    }
  }
  //cout << "objective values are the same? -" << endl;
  //for (int i = 0; i < objective_values.size() - 1; i ++)
    //assert(objective_values[i] == objective_values[i + 1]);
  write_output_file(objective_values, times, gaps, argv[1], argv[2],
    "table_COLORED.dat");

}
