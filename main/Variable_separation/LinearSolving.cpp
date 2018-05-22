using namespace std;

class LinearSolving {
  public:
	  vector<vector<int> > adjacency_list;
	  vector<vector<bool> > adjacency_matrix;
	  vector<int> f, w;
	  int N, M;
    model model_chosen;
		IloNum start_time;
	  //int initial_ub; //initial lower bound (cost of a valid tour)
	  IloEnv env;
		IloModel c_model;
		IloCplex cplex;
		/*public variables... workaround for access them in the callbacks*/
		IloNumVarArray s, e;
		IloNumArray _s, _e;
		S_cutter* s_cutter;
    vector<bool> ub_vertices;
    int ub;
    IloExpr expr_obj_fun;
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
    void startAlg (model _model_chosen);
    void endAlg (int &objective_value, double &times, double &gap);
    void Texception (int &objective_value, double &times, double &gap);
    void initial_constraints_v_infection (int initial_v);
    void add_initial_constraints_neighbors();
};


int PCI_solver::setModelProblemTopologicalModel () {
  s = IloNumVarArray(env, N, 0.0, 1.0, ILOFLOAT);
  e = IloNumVarArray(env, M, 0.0, 1.0, ILOFLOAT);
  _s = IloNumArray(env, N);
  _e = IloNumArray(env, M);
  expr_obj_fun = IloExpr(env);
  for (int v = 0; v < N; v ++) {
    sprintf(var_name, "x(%d)", v);
    s[v].setName(var_name);
    expr_obj_fun += s[v];
  }
  objective_function = IloAdd(c_model, IloMinimize(env, expr_obj_fun));
}
