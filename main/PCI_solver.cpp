
#include "PCI_solver.h"


//--------------- Constructor
PCI_solver::PCI_solver(IloEnv& _env, vector<vector<int>>& _adjacency_list,
	vector<int>& _f, vector<int>& _w) {
  // initializes cplex enviroment and model
  env = _env;
  c_model = IloModel(_env);
  cplex = IloCplex(c_model);
  //rlModel = false;
  adjacency_list = _adjacency_list;
  f = _f;
  w = _w;
	s_cutter = new S_cutter(adjacency_list, f, w);
  separation = new Separation(adjacency_list, f, w);

  N = f.size();
  ub_vertices = vector<bool>(N, true);
  ub = 0;
  for (auto weight: w) ub += weight;
  lazycall_counter = 0;
  usercall_counter = 0;
	constraints_counter = 0;
}

PCI_solver::~PCI_solver () {
  cout << "pci destructor" << endl;
  delete s_cutter;
  delete separation;
}

void PCI_solver::reset_model() {
  c_model = IloModel(env);
  cplex = IloCplex(c_model);
}

IloExpr PCI_solver::create_expression (vector<bool> selected_var) {
  IloExpr expr(env);
  for (int i = 0; i < N; i ++) {
    if (selected_var[i]) {
      expr += x[i];
		}
  }
  return expr;
}

void PCI_solver::setModelVariables(bool float_option) {
  if (float_option) {
    x = IloNumVarArray(env, N, 0.0, 1.0, ILOFLOAT);// ILOFLOAT 
  } else {
    x = IloNumVarArray(env, N, 0.0, 1.0, ILOINT);// ILOFLOAT 
  }
  _x = IloNumArray(env, N);
}

// sets model's initial constraints and objective function
void PCI_solver::setModelProblem (bool upper_bound) {
  //x = IloNumVarArray(env);
  //constraints = IloRangeArray(env);
  
  constraints = IloRangeArray(env);
  char var_name[32];

  // Objective function
  expr_obj_fun = IloExpr(env);
  for (int i = 0; i < N; i ++) {
    sprintf(var_name, "x(%d)", i);
    x[i].setName(var_name);
    expr_obj_fun += w[i] * x[i];
  }
  objective_function = IloAdd(c_model, IloMinimize(env, expr_obj_fun));

  if (upper_bound) {
    int initial_upper_bound = 
      min(s_cutter->initial_heuristic(), s_cutter->greedy_heuristic());
    IloRange ctrnt = expr_obj_fun <= initial_upper_bound;
    ctrnt.setName("heuristicUpperBoundFound");
    constraints.add(ctrnt);
  } 

  // first constraint
  vector<bool> v_aux(N, true);
  IloExpr expr_first_const = create_expression(v_aux);
  IloRange ctrnt = expr_first_const >= (*min_element(f.begin(), f.end()));
  ctrnt.setName("#1_constraint_fmin");
  constraints.add(ctrnt);

	// adds n constraints using s, and initial v = 0
  //initial_constraints_v_infection(0);
  add_initial_constraints_neighbors();
  c_model.add(constraints);
  cplex.exportModel ("model.lp");
	// to see the model we created
  //*/
}

// solving the problem after initial setup
void PCI_solver::solveProblem () {
  try {
    cplex.solve();
    cplex.getValues(_x, x);
      //cout << _x << endl;
  } catch (IloException& ex) {
    cout << "solveProblem:" << ex << endl;
  }
}

//93 sets all initial parameters
void PCI_solver::startAlg(model _model_chosen) {
	model_chosen = _model_chosen;

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
  cout << "Normal ending" << endl;
  cout << "Final obj value: " << cplex.getObjValue() << endl;
	cout << "Final gap: " << cplex.getMIPRelativeGap() << endl;
	cout << "Time: " << cplex.getTime() << endl;
  objective_value += (int) round(cplex.getObjValue());
  gap = max(min(cplex.getMIPRelativeGap(), (double) (ub - bestSol) / ub), gap);
  //gap = max(min(cplex.getMIPRelativeGap(),(ub - bestSol) / ub), gap);
  times += cplex.getTime();
  /*
  cplex.exportModel ("lpex1.lp");
  //*/
    // Get the current x solution
  /*IloNumArray xSol(env);
  cplex.getValues(xSol, x);
  vector<double> infected(N);
  cout << N << endl;
  for (int i = 0; i < N; i ++) {
    infected[i] = xSol[i];
    cout << i << ",   adj list (";
    for (auto u: adjacency_list[i]) {
      cout << u << ", ";
    }

    cout << "),   f " << f[i] << "    weight " << infected[i] << endl;
    
  }
  cout << endl;
  if (separation->finds_constraints(infected)) {
    for (int i = 0; i < (separation->constraints_rhs_res).size(); i ++) {
      for (int j = 0; j < (separation->constraints_lhs_res[i]).size(); j ++){
        cout << " v[" << separation->constraints_lhs_res[i][j] << "] "<< infected[separation->constraints_lhs_res[i][j]] << " ";
      }
      cout << " > " << separation->constraints_rhs_res[i] << endl;
    }
  }
    else {
    cout << "No constraints" << endl;
    } //*/
}


void PCI_solver::add_constraints(vector<vector<int>>& lhs, vector<int>& rhs) {
  IloRangeArray c = IloRangeArray(env);
  for (int i = 0; i < rhs.size(); i ++) {
      constraints_counter ++;
      IloExpr cutLhs(env);
      for (int j = 0; j < lhs[i].size(); j ++) {
        int k = lhs[i][j];
        cutLhs += x[k];
        //cout << k << endl;
      }
      //print_matrix(s_cutter->constraints_lhs_res, "added constraint: ");
      IloRange ctrnt = cutLhs >= rhs[i];
      //cout << ctrnt << endl;
      //ctrnt.setName("#!!!!!!!!!!!!!!!!!!!!ADDED constraint");
      c.add(ctrnt);

    }
    c_model.add(c);
}


bool PCI_solver::recursive_relax(vector<vector<int>>& lhs, vector<int>& rhs) {
  //cout << "Relaxed obj value: " << cplex.getObjValue() << endl;
  //cout << "Final gap: " << cplex.getMIPRelativeGap() << endl;
  //cout << "Time: " << cplex.getTime() << endl;
  /*
  cplex.exportModel ("lpex1.lp");
  objective_value += (int) round(cplex.getObjValue());
  gap = max(cplex.getMIPRelativeGap(), gap);
  times += cplex.getTime();
  //*/
    // Get the current x solution
  IloNumArray xSol(env);
  cplex.getValues(xSol, x);
  vector<double> infected(N);
  for (int i = 0; i < N; i ++) {
    infected[i] = xSol[i];
  }
  if (separation->finds_constraints(infected)) {
    for (int i = 0; i < (separation->constraints_rhs_res).size(); i ++) {
      lhs.push_back(separation->constraints_lhs_res[i]);
      rhs.push_back(separation->constraints_rhs_res[i]);
    }
    constraints = IloRangeArray(env);
    add_constraints(separation->constraints_lhs_res, separation->constraints_rhs_res);
    return true;
  }
  return false;
}

struct Pair {
  int v; 
  double val;
} p;

bool compare(Pair p1, Pair p2) {
  return p1.val > p2.val;
}

void relax_infected_and_selection_order(vector<bool>& infected, vector<int>& order,
  IloNumArray& _x) {
  vector<Pair> order_values;
  for (int i = 0; i < infected.size(); i ++) {
    int val = floor(_x[i] + 0.01);
    assert(val <= 1 && val >= 0);
    infected[i] = (bool) val;
    // if x[i] == 0 or fraction, can be selected
    if (0.0001 < _x[i] &&  _x[i] < 0.9999) {
      p.v = i; p.val = _x[i];
      order_values.push_back(p);
    }
  }
  sort(order_values.begin(), order_values.end(), compare);
  for (auto& p: order_values) {
    order.push_back(p.v);
  }
}

void integer_infected(vector<bool>& infected, IloNumArray& _x) {
  for (int i = 0; i < infected.size(); i ++) {
    int val = ceil(_x[i] - 0.01);
    assert(val <= 1 && val >= 0);
    infected[i] = (bool) val;
    //assert(infected[i] <= 1 && infected[i] >= 0);
  }
}

double prev_solution = 0;

bool PCI_solver::simulation_S_cutter(vector<vector<int>>& lhs, vector<int>& rhs, float considering_weight) {
  //cout << "--------__RootRelaxation-------------" << endl;
  //cout << "Relaxed obj value: " << cplex.getObjValue() << endl;
  if (cplex.getObjValue() == prev_solution) {
    return false;
  }
  prev_solution = cplex.getObjValue();
  vector<bool> infected(N);
  vector<int> order;
  bool found_new_constraints;
  
  if (considering_weight) {
    vector<double> w(N);
    for (int i = 0; i < N; i ++) {
      w[i] = _x[i];
    }
    relax_infected_and_selection_order(infected, order, _x);
    found_new_constraints = s_cutter->find_S_smaller_new_constraints(infected, model_chosen, order, w);
  } else {
    integer_infected(infected, _x);
    found_new_constraints = s_cutter->find_S_smaller_new_constraints(infected, model_chosen);
  }
  
  if (found_new_constraints) {
      // if found better upper bounds, updates ub and ub_vertices
    for (int i = 0; i < (s_cutter->constraints_rhs_res).size(); i ++) {

      lhs.push_back(s_cutter->constraints_lhs_res[i]);
      rhs.push_back(s_cutter->constraints_rhs_res[i]);
    }
    constraints = IloRangeArray(env);
    add_constraints(s_cutter->constraints_lhs_res, s_cutter->constraints_rhs_res);
    return true;
  }
  return false;
}


// exception
void PCI_solver::Texception(int &objective_value, double &times, double &gap) {
  //double elapsedTime = timer.getTime();
  //printf("%s ; %.2f ; %.2f ; %.2f ; %.2f ; %d ; %d ; %d ; %d ; %d ; %d ; %d ;
  //%d ; %d ; %d\n", filename.data(), cplex.getBestObjValue(), 100 *
  //(1 - root / ub), 100 * (1 - cplex.getBestObjValue() / ub), elapsedTime,
  //(int) cplex.getNnodes(), secCounter, precCounter, lifoCounter, sCapCounter,
  //iCapCounter, pCapCounter, rCapCounter, userCalls, lazyCalls);
  cout << "cutted off solution" << endl;
  cout << "Final obj value: " << ub << endl;
  cout << "Final gap: (cutted off result) " << cplex.getMIPRelativeGap() << 
    " or " <<  (ub - bestSol) / ub  <<endl;
  cout << "Time: " << cplex.getTime() << endl;
  //cplex.exportModel ("lpex1.lp");
  objective_value += ub;
  gap = max(min(cplex.getMIPRelativeGap(), (double) (ub - bestSol) / ub), gap);
  times += cplex.getTime();
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

// addes a constraint to initial model for each pair of vertices neighbors and
//that have trashold equal to num of neighbors
void PCI_solver::add_initial_constraints_neighbors() {
  for (int v = 0; v < N; v ++) {
    for (auto u: adjacency_list[v]) {
      if (v < u && f[v] == adjacency_list[v].size()
       && f[u] == adjacency_list[u].size()) {
        IloRange new_constraint = x[u] + x[v] >= 1;
        new_constraint.setName("#Neighbors_maxed");
        constraints.add(new_constraint);
        //add(new_constraint).end();
      }
    }
  }
}

// lazycallback adds new constraints depending on the current integer solution
// founds in the program
ILOLAZYCONSTRAINTCALLBACK1(lazyCallback, PCI_solver&, obj){
  IloInt i;
  IloEnv masterEnv = getEnv();
  int N = obj.N; //obj.x.getSize();
  // Get the current x solution
  obj.bestSol = getBestObjValue();
  cout << "---------------LazyCallback-----------------" << endl;
  cout << "Integer obj value: " << getObjValue() << " best " << obj.bestSol << " upb " << obj.ub << endl;
  cout << "gap: " << getMIPRelativeGap() << " my gap " << double (obj.ub - obj.bestSol) / obj.ub << endl;
  cout << "time: " << obj.cplex.getTime() << endl;
  IloNumArray xSol(masterEnv);
	vector<bool> infected(N, false);
  getValues(xSol, obj.x);
	for (int i = 0; i < N; i ++) {
		infected[i] = ((int) round(xSol[i])) == 1;
	}
	bool found_constraints = false;
	S_cutter* s_cutter = obj.s_cutter;
  int count_infected = 0;

  #ifdef FILE_S_CUTTER_INFO
    ofstream outfile("out_lazyconstraint.txt", ios::app);
    outfile << obj.lazycall_counter ++ << " LAZY callback:\n";
  #endif
	if (s_cutter->finds_constraints(infected, obj.model_chosen)) {
    // if found better upper bounds, updates ub and ub_vertices
          //self.tolerances.uppercutoff = min(s_cu)
		for (int i = 0; i < (s_cutter->constraints_rhs_res).size(); i ++) {
			obj.constraints_counter ++;
			IloExpr cutLhs(masterEnv);
			for (int j = 0; j < (s_cutter->constraints_lhs_res[i]).size(); j ++){
				int k = (s_cutter->constraints_lhs_res)[i][j];
				cutLhs += obj.x[k];
			}
      //print_matrix(s_cutter->constraints_lhs_res, "added constraint: ");

      #ifdef FILE_S_CUTTER_INFO
        outfile << "\nconstraint: " << cutLhs << " >= "
          << (s_cutter->constraints_rhs_res)[i];
			#endif

      add(cutLhs >= (s_cutter->constraints_rhs_res)[i]).end();
		}
	}
	else {
    // Asserts it is a valid solution found, if did not added constraints
    cout << "Solution FOUND!" << endl;
		vector<int> new_f = s_cutter->infect_graph(infected);
		for (int i = 0; i < N; i ++)
			assert(new_f[i] < 0);
	}
  #ifdef UPPERBOUND_CUT
    if (obj.ub > s_cutter->ub_val) {
      IloNumArray ub_vals(obj.env, N);
      for (int i = 0; i < N; i ++) {
        ub_vals[i] = s_cutter->current_ub[i];
      }
      obj.cplex.addMIPStart(obj.x, ub_vals, IloCplex::MIPStartCheckFeas);

      obj.ub = s_cutter->ub_val;
      //obj.ub_vertices = s_cutter->current_ub;
      //obj.cplex.setParam(IloCplex::CutUp, obj.ub - 1);
      //add(obj.expr_obj_fun <= obj.ub - 1);
      cout << "New UpperBound " << obj.ub << endl;
    }
  #endif


  #ifdef FILE_S_CUTTER_INFO
  outfile << "\n\n";
    outfile << "Objective Value: " << getObjValue() <<"  Gap: "
      << getMIPRelativeGap() << "   getBestObjValue:  "<< getBestObjValue()<< "   UpperBound Found: " << s_cutter->ub_val <<
      "  Global UpperBound: " << obj.ub << "  IloCplex UpperBound: " << obj.cplex.getParam(IloCplex::CutUp) <<  " \n";
    outfile.close();
  #endif
}



ILOUSERCUTCALLBACK1(userCallback, PCI_solver&, obj) {
  //cout << "--------UserCut-------------" << endl;
  //cout << "Obj value: " << getObjValue() << endl;
  IloEnv masterEnv = getEnv();
  int N = obj.N; //obj.x.getSize();
  // Get the current x solution
  IloNumArray xSol(masterEnv);
  getValues(xSol, obj.x);
  vector<bool> infected(N);
  vector<int> order;
  vector<double> w(N);
  for (int i = 0; i < N; i ++) {
    w[i] = xSol[i];
  }
  #ifdef FILE_S_CUTTER_INFO
    ofstream outfile("out_lazyconstraint.txt", ios::app);
    outfile << obj.usercall_counter ++ << " USER callback:\n";
  #endif
  
  vector<Pair> order_values;
  for (int i = 0; i < infected.size(); i ++) {
    int val = floor(xSol[i] + 0.01);
    assert(val <= 1 && val >= 0);
    infected[i] = (bool) val;
    if (0.0001 < xSol[i] &&  xSol[i] < 0.9999) {
      p.v = i; p.val = xSol[i];
      order_values.push_back(p);
    }
  }
  sort(order_values.begin(), order_values.end(), compare);
  for (auto& p: order_values) {
    order.push_back(p.v);
  }
  S_cutter* s_cutter = obj.s_cutter;
  bool found_new_constraints = s_cutter->find_S_smaller_new_constraints(infected, S_SMALLER_NEW, order, w);
  
  if (found_new_constraints) {
    // if found better upper bounds, updates ub and ub_vertices
    for (int i = 0; i < (s_cutter->constraints_rhs_res).size(); i ++) {
      IloExpr cutLhs(masterEnv);
      for (int j = 0; j < (s_cutter->constraints_lhs_res[i]).size(); j ++){
        int k = (s_cutter->constraints_lhs_res)[i][j];
        cutLhs += obj.x[k];
      }
      #ifdef FILE_S_CUTTER_INFO
        outfile << "\nconstraint: " << cutLhs << " >= "
          << (s_cutter->constraints_rhs_res)[i];
      #endif
      //print_matrix(s_cutter->constraints_lhs_res, "added constraint: ");
      add(cutLhs >= (s_cutter->constraints_rhs_res)[i]).end();

    }
  } 
  #ifdef UPPERBOUND_CUT
  if (obj.ub > s_cutter->ub_val) {
    IloNumArray ub_vals(obj.env, N);
    for (int i = 0; i < N; i ++) {
      ub_vals[i] = s_cutter->current_ub[i];
    }
    obj.cplex.addMIPStart(obj.x, ub_vals, IloCplex::MIPStartCheckFeas);

    obj.ub = s_cutter->ub_val;
    /*obj.ub_vertices = s_cutter->current_ub;
    obj.cplex.setParam(IloCplex::CutUp, obj.ub - 1);
    //add(obj.expr_obj_fun <= obj.ub - 1);
    //*/
    cout << "New UpperBound " << obj.ub << endl;
  }
  #endif

   #ifdef FILE_S_CUTTER_INFO
  outfile << "\n\n";
    outfile << "Objective Value: " << getObjValue() <<"  Gap: "
      << getMIPRelativeGap() << "   getBestObjValue:  "<< getBestObjValue()<< "   UpperBound Found: " << s_cutter->ub_val <<
      "  Global UpperBound: " << obj.ub << "  IloCplex UpperBound: " << obj.cplex.getParam(IloCplex::CutUp) <<  " \n";
    outfile.close();
  #endif


} 

// set cplex parameters
void PCI_solver::setCplexSettings(int timelimit) {
  #ifdef USERCUT
    cplex.use(userCallback(env, *this));
  #endif
  cplex.use(lazyCallback(env, *this));
  cout << "TIMELIMIT--------------------" << timelimit << endl;
  cplex.setParam(IloCplex::TiLim, timelimit);
  cplex.setParam(IloCplex::RandomSeed, 31415);
  cplex.setOut(env.getNullStream());
  // mip enphasis 3: CPX_MIPEMPHASIS_BESTBOUND  Emphasize moving best bound
  cplex.setParam(IloCplex::Param::Emphasis::MIP, 1);
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
    //cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, 3);
    // Strong branching
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
