/**
  PCI_solver.h
  Given an instance and a problem, using cplex api, solves the instance for and initial
  
**/
#ifndef PCISOLVER
#define PCISOLVER

#include <ilcplex/ilocplex.h>
#include <vector>
#include "Separation.h"
#include "S_cutter.h"

using namespace std;
 
class PCI_solver {
  public:
		//bool rlModel;
	  vector<vector<int>> adjacency_list;
	  vector<int> f, w;
	  int N, constraints_counter, lazycall_counter, usercall_counter, bestSol;
    model model_chosen;
		IloNum start_time;
	  //int initial_ub; //initial lower bound (cost of a valid tour)
	  IloEnv env;
		IloModel c_model;
		IloCplex cplex;
		/*public variables... workaround for access them in the callbacks*/
		IloNumVarArray x;
		IloNumArray _x;
		S_cutter* s_cutter;
    Separation* separation;

    vector<bool> ub_vertices;
    vector<bool> active;
    int ub;
    IloExpr expr_obj_fun;
  	IloObjective objective_function;
  	IloRangeArray constraints;
    // ----------------------methods
    PCI_solver (IloEnv& _env, vector<vector<int>>& _adjacency_list, vector<int>&
      _f, vector<int>& _w);
    ~PCI_solver ();
    void reset_model();
    void setModelVariables(bool float_option);
    void setModelProblem (bool upper_bound = false);
    void reduction_neighbors();
  	void solveProblem ();
    IloExpr create_expression (vector<bool>);
    void setCplexSettings (int timelimit);
    void startAlg (model _model_chosen);
    void add_constraints(vector<vector<int>>& lhs, vector<int>& rhs);
    bool recursive_relax(vector<vector<int>>& lhs, vector<int>& rhs);
    bool simulation_S_cutter(vector<vector<int>>& lhs, vector<int>& rhs, 
      float considering_weight);
    void endAlg (int &objective_value, double &times, double &gap);
    void Texception (int &objective_value, double &times, double &gap);
    void initial_constraints_v_infection (int initial_v);
    void add_initial_constraints_neighbors();
};

#endif