#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <ilcplex/ilocplex.h>
#include "S_cutter.h"

using namespace std;

class PCI_solver {
  public:
		//bool rlModel;
	  vector<vector<int> > adjacency_list;
	  vector<vector<bool> > adjacency_matrix;
	  vector<int> f, w;
	  int N, M, option, constraints_counter, lazycall_counter;
		IloNum start_time;
	  //int initial_ub; //initial lower bound (cost of a valid tour)
	  IloEnv env;
		IloModel model;
		IloCplex cplex;
		/*public variables... workaround for access them in the callbacks*/
		IloNumVarArray x;
		IloNumArray _x;
		S_cutter* s_cutter;
	  //IloNum root;
	  //IloNum ub;
  	IloObjective objective_function;
  	IloRangeArray constraints;
    // ----------------------methods
    PCI_solver (IloEnv _env, vector<vector<int> > _adjacency_list, vector<int>
      _f, vector<int> _w);
    ~PCI_solver ();
    virtual int setModelProblem ();
  	virtual int solveProblem ();
    virtual IloExpr create_expression (vector<bool>);
    void setCplexSettings (int timelimit);
    void startAlg (int _option);
    void endAlg (int &objective_value, double &times, double &gap);
    void Texception ();
    void initial_constraints_v_infection (int initial_v);
};

//--------------- Constructor
PCI_solver::PCI_solver(IloEnv _env, vector<vector<int> > _adjacency_list,
	vector<int> _f, vector<int> _w) {
  // initializes cplex enviroment and model
  env = _env;
  model = IloModel(_env);
  cplex = IloCplex(model);
  //rlModel = false;
  adjacency_list = _adjacency_list;
  f = _f;
  w = _w;
	s_cutter = new S_cutter(adjacency_list, f);
  N = f.size();
  lazycall_counter = 0;
	constraints_counter = 0;
  adjacency_matrix = vector<vector<bool> >(N, vector<bool>(N, false));
  for (int v = 0; v < N; v ++) {
    for (int i = 0; i < adjacency_list[v].size(); i ++) {
      int u = adjacency_list[v][i];
      adjacency_matrix[v][u] = true;
      adjacency_matrix[u][v] = true;
    }
  }
}

PCI_solver::~PCI_solver () {}

IloExpr PCI_solver::create_expression (vector<bool> selected_var) {
  IloExpr expr(env);
  for (int i = 0; i < N; i ++) {
    if (selected_var[i]) {
      expr += x[i];
		}
  }
  return expr;
}

// sets model's initial constraints and objective function
int PCI_solver::setModelProblem () {
  //x = IloNumVarArray(env);
  //constraints = IloRangeArray(env);
  x = IloNumVarArray(env, N, 0.0, 1.0, ILOINT);
  _x = IloNumArray(env, N);

  constraints = IloRangeArray(env);
  char var_name[32];

  // Objective function
  IloExpr expr_obj_fun(env);
  for (int i = 0; i < N; i ++) {
    sprintf(var_name, "x(%d)", i);
    x[i].setName(var_name);
    expr_obj_fun += w[i] * x[i];
  }
  objective_function = IloAdd(model, IloMinimize(env, expr_obj_fun));
  expr_obj_fun.end();

   // first constraint
  vector<bool> v_aux(N, true);
  IloExpr expr_first_const = create_expression(v_aux);
  IloRange ctrnt = expr_first_const >= (*min_element(f.begin(), f.end()));
  ctrnt.setName("#1_constraint_fmin");
  constraints.add(ctrnt);

	// adds n constraints using s, and initial v = 0
  /*initial_constraints_v_infection(0);
  model.add(constraints);
	// to see the model we created
  //*/
}

// solving the problem after initial setup
int PCI_solver::solveProblem () {
  try {
    cplex.solve();
    cplex.getValues(_x, x);
  } catch (IloException& ex) {
    cout << "solveProblem:" << ex << endl;
  }
  return 0;
}

// sets all initial parameters
void PCI_solver::startAlg(int _option) {
	option = _option;

  //nsecs = nprecs = npredCuts = nsuccCuts = nprecBal = 0;
  //root = 0.0;
  //rlModel = true;
  //ub = initial_ub;
  //#ifdef output_cuts
  //  file.open("cuts.dat");
  //#endif
  //timer.start();
	start_time = cplex.getTime();
}

// what we need done by the end of the algorithm
void PCI_solver::endAlg(int &objective_value, double &times, double &gap) {
  cout << "Final obj value: " << cplex.getObjValue() << endl;
	cout << "Final gap: " << cplex.getMIPRelativeGap() << endl;
	cout << "Time: " << cplex.getTime() << endl;
  //cplex.exportModel ("lpex1.lp");
  objective_value += (int) round(cplex.getObjValue());
  gap = max(cplex.getMIPRelativeGap(), gap);
  times += cplex.getTime();
}

// exception
void PCI_solver::Texception() {
  //double elapsedTime = timer.getTime();
  //printf("%s ; %.2f ; %.2f ; %.2f ; %.2f ; %d ; %d ; %d ; %d ; %d ; %d ; %d ;
  //%d ; %d ; %d\n", filename.data(), cplex.getBestObjValue(), 100 *
  //(1 - root / ub), 100 * (1 - cplex.getBestObjValue() / ub), elapsedTime,
  //(int) cplex.getNnodes(), secCounter, precCounter, lifoCounter, sCapCounter,
  //iCapCounter, pCapCounter, rCapCounter, userCalls, lazyCalls);
}

// lazycallback adds new constraints depending on the current integer solution
// founds in the program
ILOLAZYCONSTRAINTCALLBACK1(lazyCallback, PCI_solver&, obj){
  IloInt i;
  IloEnv masterEnv = getEnv();
  int N = obj.N; //obj.x.getSize();
  // Get the current x solution

  IloNumArray xSol(masterEnv);
	vector<bool> infected(obj.N, false);
  getValues(xSol, obj.x);
	for (int i = 0; i < obj.N; i ++) {
		infected[i] = ((int) round(xSol[i])) == 1;
	}
	bool found_constraints = false;
	S_cutter* s_cutter = obj.s_cutter;
  int count_infected = 0;

  ofstream outfile("out_lazyconstraint.txt", ios::app);
  outfile << obj.lazycall_counter ++ << " callback:\n";
  outfile << "Objective Value: " << getObjValue() << "\n Infected vertices: "
    << infected[0];
  for (int i = 1 ; i < obj.N; i ++) {
    outfile << ", " << infected[i];
  }


	if (s_cutter->finds_constraints(infected, obj.option)) {
          //self.tolerances.uppercutoff = min(s_cu)
		for (int i = 0; i < (s_cutter->constraints_rhs_res).size(); i ++) {
			obj.constraints_counter ++;
			IloExpr cutLhs(masterEnv);
			for (int j = 0; j < (s_cutter->constraints_lhs_res[i]).size(); j ++){
				int k = (s_cutter->constraints_lhs_res)[i][j];
				cutLhs += obj.x[k];
			}
      //print_matrix(s_cutter->constraints_lhs_res, "added constraint: ");
      outfile << "\nconstraint: " << cutLhs << " >= " << (s_cutter->constraints_rhs_res)[i];
			add(cutLhs >= (s_cutter->constraints_rhs_res)[i]).end();
		}
	}
	else {
    // Asserts it is a valid solution found, if did not added constraints
		vector<int> new_f = s_cutter->infect_graph(infected);
		for (int i = 0; i < N; i ++)
			assert(new_f[i] < 0);
	}
  outfile << "\n";
  outfile.close();
}


// added initial constraints v, so starts with a more robust model
void PCI_solver::initial_constraints_v_infection(int initial_v) {
  vector<bool> infected (N, false);//, used_vertices (N, false);
  int used_count = 0;
  infected[initial_v] = true;
  //used_vertices[initial_v] = true;
  //S_cutter* s_cutter = s_cutter;
  // adds at least N constraints at the initial model
  while (used_count ++ < N) {
    if (s_cutter->finds_s_smaller_constraints(infected, S_SMALLER)) {
      // will generate a new infected vector, for the next iteration
      //fill(infected.begin(), infected.end(), false);
      //for (int i = 0; i < (s_cutter->rhs).size(); i ++) {
      IloExpr cutLhs(env);
      int k = 0;
      initial_v = 0;
      for (int j = 0; j < (s_cutter->constraints_lhs_res[0]).size(); j ++){
        k = (s_cutter->constraints_lhs_res)[0][j];
        cutLhs += x[k];
        if (! infected[k])
          initial_v = k;
      }
      if (infected[initial_v]) {
        for (int i = 0; i < N; i ++)
          if (! infected[i])
            initial_v = i;
      }
      infected[initial_v] = true;
      //used_count ++;
      //copy(used_vertices.begin(), used_vertices.end(), infected.begin());
      IloRange ctrnt = cutLhs >= (s_cutter->constraints_rhs_res)[0];
      ctrnt.setName("#S_initial_constraint");
      constraints.add(ctrnt);
      //add(cutLhs >= (s_cutter->rhs)[i]).end();
      //}
    }
  }
}

// set cplex parameters
void PCI_solver::setCplexSettings(int timelimit) {
  cplex.use(lazyCallback(env, *this));
  cplex.setParam(IloCplex::TiLim, timelimit);
  cplex.setParam(IloCplex::RandomSeed, 31415);
  cplex.setOut(env.getNullStream());
  /*cplex.setParam(IloCplex::CutUp, initial_ub); //Sets the upper cutoff tolerance.

    if (vlDisp != 0) {
      cplex.setParam(IloCplex::MIPDisplay, 2);
      cplex.setParam(IloCplex::MIPInterval, vlDisp);
    } else
      cplex.setOut(Env.getNullStream());
    cplex.setWarning(Env.getNullStream());
    //*/
    /*cplex.setParam(IloCplex::RootAlg, alg);
    cplex.setParam(IloCplex::NodeAlg, alg);
    cplex.setParam(IloCplex::Threads, numThreads);//*/

    /*if (vlEmph != 0)
      cplex.setParam(IloCplex::MIPEmphasis, vlEmph);
    if (vlGap != 0.0) {
      cplex.setParam(IloCplex::EpAGap, 1 - vlGap);
      cplex.setParam(IloCplex::EpInt, vlGap);
      cplex.setParam(IloCplex::ObjDif, 1 - vlGap);
      //cplex.setParam(IloCplex::EpOpt,vlGap);
      //cplex.setParam(IloCplex::EpGap,vlGap);
    }

    if (mem != 0)
      cplex.setParam(IloCplex::WorkMem, mem);

    cplex.setParam(IloCplex::MemoryEmphasis, 1);

    //cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, 3); // Strong branching
    //cplex.setParam(IloCplex::VarSel, 3);
    //cplex.setParam(IloCplex::StrongItLim, 1);

    cplex.setParam(IloCplex::WorkDir, ".");
    cplex.setParam(IloCplex::NodeFileInd, 2);
    //TreLim
    //*//*
    cplex.setParam(IloCplex::HeurFreq, -1); // heuristic frequency
    cplex.setParam(IloCplex::RINSHeur, -1); // do not apply RINS heuristic
    cplex.setParam(IloCplex::FPHeur, -1); // do not apply feasibility pump heuristic
    cplex.setParam(IloCplex::LBHeur, 0); // do not apply local branching heuristic
    cplex.setParam(IloCplex::PreInd, 0); // do not apply presolve
    cplex.setParam(IloCplex::PreslvNd, -1); // node presolve
    cplex.setParam(IloCplex::Symmetry, 0); // symmetry breaking
    cplex.setParam(IloCplex::AggInd, 0); // do not use any aggregator
    cplex.setParam(IloCplex::BndStrenInd, 0); // no var bound strengthening
    cplex.setParam(IloCplex::CoeRedInd, 0); // no coefficient reduction
    cplex.setParam(IloCplex::DepInd, 0); // no dependency checker
    cplex.setParam(IloCplex::Reduce, CPX_PREREDUCE_NOPRIMALORDUAL); // no reductions
    cplex.setParam(IloCplex::CutPass, -1); // cutting plane passes at root node
    cplex.setParam(IloCplex::Cliques, -1);
    cplex.setParam(IloCplex::Covers, -1);
    cplex.setParam(IloCplex::DisjCuts, -1);
    cplex.setParam(IloCplex::FlowCovers, -1);
    cplex.setParam(IloCplex::FlowPaths, -1);
    cplex.setParam(IloCplex::FracCuts, -1);
    cplex.setParam(IloCplex::GUBCovers, -1);
    cplex.setParam(IloCplex::ImplBd, -1);
    cplex.setParam(IloCplex::MIRCuts, -1);
    cplex.setParam(IloCplex::MCFCuts, -1);
    cplex.setParam(IloCplex::ZeroHalfCuts, -1);
  /*/
}
