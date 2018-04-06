#ifndef TSPPDMS_H
#define TSPPDMS_H

#include <iostream>
#include <vector>
#include <map>
#include <stack>
#include <algorithm>
#include <cstdarg>
#include <ilcplex/ilocplex.h>
#include "Dinic.h"

using namespace std;

#define MAXNODES 50
#define MAX_ID 1000
#define MAXSTKS 4

typedef struct {
	int u, v; //
	int c; //original edge cost
} Arc;

class TSPPDMS {
private:
	int V, A, n; //nodes, arcs and requisition num.
	vector< pair<float, float> > nodesCoord; //only for instances given as coordinates
	vector<Arc*> arcs; //graph edges
	vector<Arc*> adj[MAXNODES]; //adjacency list for each node

	int initial_ub; //initial upper bound (cost of a valid tour)
	string filename;

	IloObjective fo;
	IloRangeArray Constraints;

	IloConversion convy;
	IloConversion convx;

	int path[MAXNODES];

	int nsecs;
	int nprecs;
	int npredCuts;
	int nsuccCuts;
	int nprecBal;

	IloNum root;
	IloNum ub;
	bool rlModel;

public:

	TSPPDMS(IloEnv e);
	virtual int readData(const char[]);
	virtual int setModelProblem();
	virtual int solveProblem();
	virtual void relaxIntVars();
	virtual void enforceIntVars();
	virtual IloNum getSolution();
	virtual void addCut(IloRange);
	virtual void addCuts();
	virtual bool findCuts(int&, int&, int&, int&);
	virtual void setCplexSettings(int, int, int, int, double, int);
	virtual void startAlg();
	virtual void endAlg();
	virtual void Texception();
	void create_graphviz_solution();
	void create_graphviz_linear_solution();
	void create_graphviz_image(int, int);


	int create_svg(const char*);
	int getNumRequisition();
	bool checkSolution();
	void setUB(IloNum);
	void modifyGraph(Dinic&, int, int, int);
	IloRange search_violated_rounded_capacity_inequalities(bool&);

	IloRange generateSEC(IloNum&, IloExpr&);
	IloRange generatePREC(IloNum&, IloExpr&);
	IloRange generatePRECI(IloNum&, IloExpr&, int);
	IloRange generatePredecessorCut(IloNum&, IloExpr&);
	IloRange generateSuccessorCut(IloNum&, IloExpr&);
	IloRange generateLIFO(int, int, int, IloNum&, IloExpr&);
	IloRange generateLIFOS(int, int, int, IloNum&, IloExpr&);
	IloRange generate_lifted_lifo(int, int, int, IloNum&, IloExpr&, vector<int>&);
	IloRange generate_lifted_lifo_loading_inside(int, int, int, IloNum&, IloExpr&);
	IloRange generateLIFOS2(int, int, int, IloNum&, IloExpr&);
	IloRange generateConflicitCapacity(int, IloNum&, IloExpr&);
	IloRange generateRoundedCapacity(IloNum&, IloExpr&);
	IloRange generateSECCapacity(int, IloNum&, IloExpr&, vector<int>&);
	IloRange generateSECCapacity2(int, int&, IloNum&, IloExpr&);

	IloEnv Env;
	IloTimer Timer;
	IloModel Model;
	IloCplex Cplex;

	/*public variables... workaround for access them in the callbacks*/
	IloNumVarArray y;
	IloNumArray _y;
	IloNumVarArray x;
	IloNumArray _x;

	map<int, int> mapVarsX;
	map<int, int> mapVarsY;

	bool is_pickup[MAXNODES]; //1-pickup, 0-delivery
	int sibling[MAXNODES]; //request (p,d) => sibling[p]=d sibling[d]=p
	vector<int> pickups;
	vector<int> deliveries;
	int Q; //stack capacity
	int K; //number of available stacks
	short int itemLength[MAXNODES];

	vector<int> nodeSubset; //defines a subset of the nodes sets
	vector<int> compSubset;

	IloRangeArray secs;
	IloRangeArray precs;
	IloRangeArray lifos;
	IloRangeArray caps;

	IloExprArray LHS;
	IloNumArray RHS;
};

#endif  /* TSPPDMS_H */
