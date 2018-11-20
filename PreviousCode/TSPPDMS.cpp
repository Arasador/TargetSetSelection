#include "TSPPDMS.h"

#define EPSILON 0.01
#define LIFTED_PRECEDENCE_SEC // LPRS
#define SUCC_INEQ           // SUCC
#define INITIAL_PRECEDENCE_S1 // IPS1
#define INITIAL_PRECEDENCE_S2 // IPS2
//#define INITIAL_PRECEDENCE_S3 // IPS3
//#define ORIGINAL_LIFO       // OL
#define LIFTED_LIFO           // LL
//#define PATH_CAPACITY       // PATHCAP
#define SEC_CAPACITY          // SECCAP
#define FIND_VIOLATED_PATH_CAPACITY_AND_LIFO_IN_FINDCUTS // PATHFC
#define FIND_VIOLATED_PATH_CAPACITY_AND_LIFO_IN_USERCALLBACK // PATHUC
//#define SEPARATE_LIFO_EXACTLY_IN_ILOUSER

#define INTEGER_SOLUTION 0
#define LINEAR_SOLUTION 1

bool lifted_lifo = true;
bool original_lifo = false;
bool sec_capacity = true;
bool path_capacity = false;
int secCounter = 0;
int precCounter = 0;
int lifoCounter = 0;
int sCapCounter = 0;
int iCapCounter = 0;
int pCapCounter = 0;
int rCapCounter = 0;
int userCalls = 0;
int lazyCalls = 0;
int fcuts = 0;
int currentNode = 0;

TSPPDMS::TSPPDMS(IloEnv e) : Env(e), Timer(e) {
  Model = IloModel(Env);
  Cplex = IloCplex(Model);
  rlModel = false;
}

void TSPPDMS::startAlg() {
  nodeSubset = vector<int>(2 * n);
  compSubset = vector<int>(2 * n);
  nodeSubset.reserve(2 * n);
  compSubset.reserve(2 * n);
  nsecs = nprecs = npredCuts = nsuccCuts = nprecBal = 0;
  //root = 0.0;
  //rlModel = true;
  ub = initial_ub;
#ifdef output_cuts
  file.open("cuts.dat");
#endif
  Timer.start();
}

void TSPPDMS::endAlg() {
  double elapsedTime = Timer.getTime();
  printf("%s ; %.2f ; %.2f ; %.2f ; %.2f ; %d ; %d ; %d ; %d ; %d ; %d ; %d ; %d ; %d ; %d\n", filename.data(), Cplex.getObjValue(), 100 * (1 - root / ub), 100 * Cplex.getMIPRelativeGap(), elapsedTime, (int) Cplex.getNnodes(), secCounter, precCounter, lifoCounter, sCapCounter, iCapCounter, pCapCounter, rCapCounter, userCalls, lazyCalls);
}

void TSPPDMS::Texception() {

  double elapsedTime = Timer.getTime();
  printf("%s ; %.2f ; %.2f ; %.2f ; %.2f ; %d ; %d ; %d ; %d ; %d ; %d ; %d ; %d ; %d ; %d\n", filename.data(), Cplex.getBestObjValue(), 100 * (1 - root / ub), 100 * (1 - Cplex.getBestObjValue() / ub), elapsedTime, (int) Cplex.getNnodes(), secCounter, precCounter, lifoCounter, sCapCounter, iCapCounter, pCapCounter, rCapCounter, userCalls, lazyCalls);

}

ILOMIPCALLBACK0(MipCallback) {
  currentNode++;
}

ILOLAZYCONSTRAINTCALLBACK1(lazyCallback, TSPPDMS&, obj) {

  short int component[MAXNODES]; //defines the connected component of each node
  short int comp; //components counter
  short int path[MAXNODES]; //defines the path sequence
  short int pIndex; //for indexing the discovered path
  short int u, v; //auxiliary
  bool flag;
  bool secFound, precFound, lifoFound, conflictFound;
  secFound = false;
  precFound = false;
  lifoFound = false;
  conflictFound = false;
  int n;
  IloNum ctValue1;
  bool visited[MAXNODES];
  bool inS[MAXNODES];
  bool precedence[MAXNODES];
  int used_stack[MAXNODES];
  bool precs[MAXNODES];
  double qs;

  IloExpr expr(getEnv());
  n = obj.getNumRequisition();

  try {
    lazyCalls++;
    getValues(obj._x, obj.x); //path
    getValues(obj._y, obj.y); //stack


    /* SUBTOUR ELIMINATION CONSTRAINTS */
    memset(component, 0, sizeof (component));
    memset(path, 0, sizeof (path));

    u = 0;
    pIndex = 0;
    flag = false;
    while (!flag) {
      for (int j = 1; j <= 2 * n + 1; j++) {
        if (u != j && !(!obj.is_pickup[u] && obj.sibling[u] == j) && j != 0 && u != 2 * n + 1
                && !(obj.is_pickup[u] && j == 2 * n + 1) && !(u == 0 && !obj.is_pickup[j])) {
          if (obj._x[ obj.mapVarsX[100 * u + j] ] > 1 - EPSILON) {
            if (obj.is_pickup[j]) {
              for (int k = 0; k < obj.K; k++) {
                if (obj._y[ obj.mapVarsY[10 * j + k] ] > 1 - EPSILON) {
                  used_stack[j] = k;
                  used_stack[obj.sibling[j]] = k;
                  k = obj.K;
                }
              }
            }
            component[u] = 1;
            u = j;
            j = 2 * n + 2;
            if (u == 2 * n + 1)
              flag = true;
            path[++pIndex] = u;
          }
        }
      }
    }

    secFound = false;
    component[u] = 1;
    comp = 1;
    for (int i = 0; i <= 2 * n + 1; i++) {
      if (component[i] == 0) {
        u = i;
        flag = false;
        comp++;
        obj.nodeSubset.clear();
        while (!flag) {
          for (int j = 0; j <= 2 * n + 1; j++) {
            if (u != j && !(!obj.is_pickup[u] && obj.sibling[u] == j) && j != 0 && u != 2 * n + 1
                    && !(obj.is_pickup[u] && j == 2 * n + 1) && !(u == 0 && !obj.is_pickup[j])) {
              if (obj._x[ obj.mapVarsX[100 * u + j] ] > 1 - EPSILON) {
                component[u] = comp;
                obj.nodeSubset.push_back(u);
                u = j;
                j = 2 * n + 2;
                if (u == i) {
                  flag = true;
                  IloRange cut;
                  ctValue1 = 0.0;
                  cut = obj.generateSEC(ctValue1, expr);
                  add(cut, IloCplex::UseCutPurge);
                  expr.clear();
                  secFound = true;
                  secCounter++;
                }
              }
            }
          }
        }
      }
    }

    /**********************************************************/
    /*************PRECEDENCE CONSTRAINTS***********************/
    /**********************************************************/
    precFound = false;
    memset(visited, false, sizeof (visited));
    memset(precedence, true, sizeof (precedence));
    obj.nodeSubset.clear();
    obj.nodeSubset.push_back(0);
    for (int i = 1; i < pIndex; i++) {
      u = path[i];
      obj.nodeSubset.push_back(u);
      v = obj.sibling[u];
      visited[u] = true;
      if (!obj.is_pickup[u] && !visited[v]) {
        IloRange cut;
        ctValue1 = 0.0;
        cut = obj.generatePREC(ctValue1, expr);
        expr.clear();
        add(cut, IloCplex::UseCutPurge);
        precFound = true;
        precedence[u] = precedence[v] = false;
        precCounter++;
      }
    }


    /************************************************************/
    /*******************LIFO CONSTRAINTS*************************/
    /************************************************************/
    lifoFound = false;
    conflictFound = false;

    u = 0;
    for (int i = 1; i < pIndex; i++) {
      v = path[i];
      if (obj.is_pickup[v] && precedence[v] && component[obj.sibling[v]] == 1) //v is a pickup node
      {
        int stk = used_stack[v];
        int dv = obj.sibling[v];
        int last = 0;
        memset(visited, false, sizeof (visited));
        memset(inS, false, sizeof (inS));

        u = v;
        for (int ii = i + 1; ii < pIndex; ii++) {
          v = path[ii];
          if (v != dv) {
            if (used_stack[v] == stk) {
              visited[v] = true;
            }
            u = v;
          } else {
            last = ii;
            ii = pIndex; //breaks the for...
          }
        }

        flag = false;
        obj.nodeSubset.clear();

#ifdef LIFTED_LIFO

        bool possible_lifo = false;
        vector<int> loaded_in_k;

        for (int j = last - 1; j > i; --j) {
          u = path[j];
          v = obj.sibling[u];
          obj.nodeSubset.push_back(u);
          inS[u] = true;
          if (obj.is_pickup[u] && visited[u] && !visited[v]) {
            loaded_in_k.push_back(u);
            j = i;
            possible_lifo = true;
          }
        }

        if (possible_lifo) {
          for (int j = 0; j < obj.nodeSubset.size(); j++) {
            u = obj.nodeSubset[j];
            v = obj.sibling[u];
            if (obj.is_pickup[u] && !inS[v] && obj.nodeSubset[j] != loaded_in_k[0]) {
              loaded_in_k.push_back(u);
            }
          }
          ctValue1 = 0.0;
          IloRange cut = obj.generate_lifted_lifo(path[i], loaded_in_k[0], stk, ctValue1, expr, loaded_in_k);
          expr.clear();
          add(cut, IloCplex::UseCutPurge);
          flag = true;
          lifoFound = true;
          lifoCounter++;
        }

        loaded_in_k.clear();
        obj.nodeSubset.clear();
        possible_lifo = false;
        memset(inS, false, sizeof (inS));

        obj.nodeSubset.push_back(path[i]);
        loaded_in_k.push_back(path[i]);

        int last_delivery = -1;
        for (int j = i + 1; j < last; ++j) {
          u = path[j];
          v = obj.sibling[u];
          obj.nodeSubset.push_back(u);
          inS[u] = true;
          if (!obj.is_pickup[u] && visited[u] && !visited[v] && !precedence[v]) {
            last_delivery = j;
            j = last;
            possible_lifo = true;
            obj.nodeSubset.pop_back();
          }
        }

        if (possible_lifo) {
          for (int j = 0; j < obj.nodeSubset.size(); j++) {
            u = obj.nodeSubset[j];
            v = obj.sibling[u];
            if (obj.is_pickup[u] && !inS[v] && obj.nodeSubset[j] != loaded_in_k[0]) {
              loaded_in_k.push_back(u);
            }
          }
          ctValue1 = 0.0;
          IloRange cut = obj.generate_lifted_lifo(obj.sibling[path[last_delivery]], path[i], stk, ctValue1, expr, loaded_in_k);
          expr.clear();
          add(cut, IloCplex::UseCutPurge);
          lifoCounter++;
          flag = true;
          lifoFound = true;
        }


        loaded_in_k.clear();
        obj.nodeSubset.clear();

#endif        


#ifdef ORIGINAL_LIFO

        for (int j = last - 1; j >= i + 1; j--) {
          u = path[j];
          v = obj.sibling[u];
          obj.nodeSubset.push_back(u);
          if (obj.is_pickup[u]) {
            if (visited[u] && !visited[v])//a pickup in S for which the delivery is outside S
            {
              //for (int kk = 0; kk < obj.K; kk++) {
              ctValue1 = 0.0;
              IloRange cut = obj.generateLIFOS(path[i], u, stk, ctValue1, expr);
              expr.clear();
              add(cut, IloCplex::UseCutPurge);
              j = 0;
              flag = true;
              lifoFound = true;
              lifoCounter++;
              //}
            }
          }
        }

        obj.nodeSubset.clear();

        obj.nodeSubset.push_back(path[i]);

        for (int j = i + 1; j < last; j++) {
          u = path[j];
          v = obj.sibling[u];
          obj.nodeSubset.push_back(u);
          if (!obj.is_pickup[u]) {
            if (visited[u] && !visited[v] && !precedence[v])//a delivery in S for which the pickup is outside S
            {
              visited[u] = false;
              obj.nodeSubset.pop_back();
              //for (int kk = 0; kk < obj.K; kk++) {
              ctValue1 = 0.0;
              IloRange cut = obj.generateLIFOS(v, path[i], stk, ctValue1, expr);
              expr.clear();
              add(cut, IloCplex::UseCutPurge);
              j = last;
              flag = true;
              lifoFound = true;
              lifoCounter++;
              //}
            }
          }
        }

        obj.nodeSubset.clear();

#endif

        if (!flag) {
          obj.nodeSubset.clear();
          for (int j = i + 1; j < last; j++)
            obj.nodeSubset.push_back(path[j]);
          /****************************************
           *******Conflict Capacity Constraints*****
           *****************************************/
          short int phi = 0;
          short int the = 0;
          short int zS = 0;
          vector<int> predecessors;
          vector<int> successors;
          memset(visited, false, sizeof (visited));
          conflictFound = false;
          qs = 0;
          for (int l = 0; l < obj.nodeSubset.size(); l++) {
            if (!obj.is_pickup[obj.nodeSubset[l]]) {
              predecessors.push_back(obj.sibling[obj.nodeSubset[l]]);
              qs -= 1; //needless
            } else {
              successors.push_back(obj.sibling[obj.nodeSubset[l]]);
              qs += 1; //needless
            }
            visited[obj.nodeSubset[l]] = true;
          }
          for (int l = 0; l < predecessors.size(); l++)
            if (!visited[predecessors[l]])
              phi += obj.itemLength[predecessors[l]];
          for (int l = 0; l < successors.size(); l++)
            if (!visited[successors[l]])
              the += -obj.itemLength[successors[l]];
          zS = phi >= the ? phi : the;
          if (zS > (obj.Q)*(obj.K - 1)) {
            ctValue1 = 0.0;
            IloRange capCut = obj.generateConflicitCapacity(path[i], ctValue1, expr);
            add(capCut, IloCplex::UseCutPurge);
            expr.clear();
            conflictFound = true;
            iCapCounter++;
          }
          obj.nodeSubset.clear();
        }
      }//if(obj.is_pickup[v])
      u = path[i];
    }

    /*********************************************
     *********Capacity Path Inequalities***********
     **********************************************/
#ifdef PATH_CAPACITY
    int count = 0;
    if (!precFound && !lifoFound) {
      short int start[MAXSTKS];
      short int end[MAXSTKS];
      short int capacities[MAXSTKS];
      short int old_capacity;

      u = 0;
      expr.clear();
      memset(capacities, 0, sizeof (capacities));
      memset(visited, false, sizeof (visited));
      for (int i = 1; i < pIndex; i++) {
        v = path[i];
        int k = used_stack[v];
        if (obj.is_pickup[v]) {
          old_capacity = capacities[k];
          capacities[k] += obj.itemLength[v];
          visited[v] = true;
          if (old_capacity == 0) {
            start[k] = i - 1;
          } else if (capacities[k] > obj.Q) {
            end[k] = i;
            //for (int kk = 0; kk < obj.K; kk++) {//generate a path inequality
            for (int ti = start[k]; ti < end[k]; ti++) {
              if (obj.is_pickup[path[ti + 1]] && used_stack[path[ti + 1]] == k) {
                expr += obj.y[ obj.mapVarsY[10 * path[ti + 1] + k] ];
                expr += obj.x[ obj.mapVarsX[100 * path[ti] + path[ti + 1]] ];
                count += 2;
              } else {
                count++;
                expr += obj.x[ obj.mapVarsX[100 * path[ti] + path[ti + 1]] ];
              }
            }
            add(expr <= count - 1, IloCplex::UseCutPurge).end();
            expr.clear();
            count = 0;
            pCapCounter++;
            //}
          }
        }
        if (!obj.is_pickup[v] && visited[obj.sibling[v]]) {
          capacities[k] -= -obj.itemLength[v];
        }
        u = v;
      }
    }
#endif

#ifdef SEC_CAPACITY
    if (!precFound && !lifoFound) {
      short int start[MAXSTKS];
      short int end[MAXSTKS];
      short int capacities[MAXSTKS];
      short int di;
      short int old_capacity;
      vector<int> loaded_in_k;

      u = 0;
      expr.clear();
      memset(capacities, 0, sizeof (capacities));
      memset(visited, false, sizeof (visited));
      for (int i = 1; i < pIndex; i++) {
        v = path[i];
        int k = used_stack[v];
        if (obj.is_pickup[v]) {
          old_capacity = capacities[k];
          capacities[k] += obj.itemLength[v];
          visited[v] = true;
          if (old_capacity == 0) //before loading v, stack k was empty...
          {
            start[k] = i - 1;
          } else if (capacities[k] > obj.Q) {
            end[k] = i;
            //for (int kk = 0; kk < obj.K; kk++) {
            obj.nodeSubset.push_back(path[start[k]]);
            for (int ti = start[k]; ti < end[k]; ti++) {
              if (obj.is_pickup[path[ti + 1]] && used_stack[path[ti + 1]] == k) {
                obj.nodeSubset.push_back(path[ti + 1]);
                loaded_in_k.push_back(path[ti + 1]);
              } else {
                obj.nodeSubset.push_back(path[ti + 1]);
              }
            }
            ctValue1 = 0.0;
            IloRange cut = obj.generateSECCapacity(k, ctValue1, expr, loaded_in_k);
            expr.clear();
            loaded_in_k.clear();
            obj.nodeSubset.clear();
            add(cut, IloCplex::UseCutPurge);
            pCapCounter++;
            //}
          }
        }
        if (!obj.is_pickup[v] && visited[obj.sibling[v]]) {
          capacities[k] -= -obj.itemLength[v];
        }
        u = v;
      }
    }
#endif

  } catch (IloException& ex) {
    cout << "lazyCallback: " << ex << endl;
  }
}

ILOUSERCUTCALLBACK1(cutCallback, TSPPDMS&, obj) {
  if (!isAfterCutLoop())
    return;

  int s, t, u, v;
  int flow, n, size;
  IloNum ctValue1, ctValue2, ctValue3;
  bool flag;
  bool incS[MAXNODES];
  bool capa[MAXNODES];
  bool secFound, precFound, capFound;

  IloExpr expr(getEnv());
  n = obj.getNumRequisition();

  try {

    userCalls++;
    getValues(obj._x, obj.x);
    getValues(obj._y, obj.y);

    for (int i = 0; i < 5; i++) {
      bool foundCut = false;
      IloRange cut = obj.search_violated_rounded_capacity_inequalities(foundCut);
      if (foundCut) {
        add(cut, IloCplex::UseCutPurge).end();
        rCapCounter++;
        i=5;
      }
    }

    Dinic suppGraph(2 * n + 2);
    for (int i = 1; i <= 2 * n; i++) {
      for (int j = 1; j <= 2 * n; j++) {
        if (i != j && (obj.is_pickup[i] || obj.sibling[i] != j)) {
          if (obj._x[ obj.mapVarsX[100 * i + j] ] > 0.0001){
            suppGraph.addArc(i, j, (int) (10000 * obj._x[ obj.mapVarsX[100 * i + j] ] + 0.5));
          }
        }
      }
    }

    for (int i = 0; i < obj.pickups.size(); i++) {
      if (obj._x[ obj.mapVarsX[100 * 0 + obj.pickups[i]] ] > 0.0001) {
        suppGraph.addArc(0, obj.pickups[i], (int) (10000 * obj._x[ obj.mapVarsX[100 * 0 + obj.pickups[i]] ] + 0.5));
      }
    }
    for (int i = 0; i < obj.deliveries.size(); i++) {
      if (obj._x[ obj.mapVarsX[100 * obj.deliveries[i]+(2 * n + 1)] ] > 0.0001) {
        suppGraph.addArc(obj.deliveries[i], 2 * n + 1, (int) (10000 * obj._x[ obj.mapVarsX[100 * obj.deliveries[i]+(2 * n + 1)] ] + 0.5));
      } 
    }
    secFound = false;
    precFound = false;
    capFound = false;
    /****************************************************************/
    /*************SUBTOUR ELIMINATION CONSTRAINTS********************/
    /****************************************************************/
    //solves the maxflow from depot 0 to each of the pickups
    s = 0;
    memset(incS, false, sizeof (incS));
    for (int i = 0; i < obj.pickups.size(); i++) {
      t = obj.pickups[i];
      if (!incS[t]) {
        Dinic graph = suppGraph;
        flow = graph.maxFlow(s, t);
        if (flow < 9990) {
          graph.runMinCut();
          obj.nodeSubset.clear();
          for (int ii = 0; ii <= 2 * n + 1; ii++)
            if (!graph.isInMinCut(ii))
              obj.nodeSubset.push_back(ii);

          if (obj.nodeSubset.size() > 1) {
            IloRange cut;
            ctValue1 = 0.0;
            cut = obj.generateSEC(ctValue1, expr);
            expr.clear();
            if(ctValue1 > obj.nodeSubset.size() - 1 + EPSILON){
              add(cut, IloCplex::UseCutPurge).end();
              secFound = true;
              secCounter++;            
              for (int ii = 0; ii < obj.nodeSubset.size(); ii++) {
                incS[obj.nodeSubset[ii]] = true;
              }
            }
          }
        } else if (!capFound) {
          graph.runMinCut();
          obj.nodeSubset.clear();
          for (int ii = 1; ii < 2 * n + 1; ii++)
            if (graph.isInMinCut(ii))
              obj.nodeSubset.push_back(ii);

          double qS = 0;
          flag = false;
          for (int ii = 0; ii < obj.nodeSubset.size(); ii++) {
            qS += obj.itemLength[obj.nodeSubset[ii]];
            if (obj.nodeSubset[ii] == 0 || obj.nodeSubset[ii] == 2 * n + 1)
              flag = true;
          }
          if (qS < 0) {
            qS = -qS;
          }
          if (qS > obj.Q && !flag) {
            ctValue1 = 0.0;
            IloRange cut = obj.generateRoundedCapacity(ctValue1, expr);
            if (ctValue1 + EPSILON < ceil(qS / (obj.K * obj.Q))) {
              add(cut, IloCplex::UseCutPurge).end();
              capFound = true;
              rCapCounter++;
            }
            expr.clear();
          }
        }//else*/
      }
    }
    //solves the maxflow from each of the deliveires to the depot 2n+1
    t = 2 * n + 1;
    memset(incS, false, sizeof (incS));
    capFound = false;
    for (int i = 0; i < obj.deliveries.size(); i++) {
      s = obj.deliveries[i];
      if (!incS[s]) {
        Dinic graph = suppGraph;
        flow = graph.maxFlow(s, t);
        if (flow < 9990) {
          graph.runMinCut();
          obj.nodeSubset.clear();
          for (int ii = 0; ii <= 2 * n + 1; ii++)
            if (graph.isInMinCut(ii))
              obj.nodeSubset.push_back(ii);
          if (obj.nodeSubset.size() > 1) {
            IloRange cut;
            ctValue1 = 0.0;
            cut = obj.generateSEC(ctValue1, expr);
            expr.clear();
            if(ctValue1 > obj.nodeSubset.size() - 1 + EPSILON){
              add(cut, IloCplex::UseCutPurge).end();
              secFound = true;
              secCounter++;            
              for (int ii = 0; ii < obj.nodeSubset.size(); ii++){
                incS[obj.nodeSubset[ii]] = true;
              }
            }
          }
        } else if (!capFound) {
          graph.runMinCut();
          obj.nodeSubset.clear();
          for (int ii = 1; ii < 2 * n + 1; ii++)
            if (graph.isInMinCut(ii))
              obj.nodeSubset.push_back(ii);
          double qS = 0;
          flag = false;
          for (int ii = 0; ii < obj.nodeSubset.size(); ii++) {
            qS += obj.itemLength[obj.nodeSubset[ii]];
            if (obj.nodeSubset[ii] == 0 || obj.nodeSubset[ii] == 2 * n + 1) {
              flag = true;
            }
          }
          if (qS < 0)
            qS = -qS;
          if (qS > obj.Q && !flag) {
            ctValue1 = 0.0;
            IloRange cut = obj.generateRoundedCapacity(ctValue1, expr);
            if (ctValue1 + EPSILON < ceil(qS / (obj.K * obj.Q))) {
              add(cut, IloCplex::UseCutPurge).end();
              capFound = true;
              rCapCounter++;
            }
            expr.clear();
          }

        }
      }
    }

    /****************************************************************/
    /*************PRECEDENCE CONSTRAINTS*****************************/
    /****************************************************************/
    for (int i = 0; i < n; i++) {
      u = obj.pickups[i];
      v = obj.sibling[u];
      s = 0;
      t = 2 * n + 1;
      Dinic graph = suppGraph;
      graph.addArc(s, v, 20000);
      graph.addArc(u, t, 20000);
      flow = graph.maxFlow(s, t);
      if (flow < 19990) {
        graph.runMinCut();
        obj.nodeSubset.clear();
        for (int ii = 0; ii <= 2 * n + 1; ii++) {
          if (graph.isInMinCut(ii)) {
            obj.nodeSubset.push_back(ii);
          }
        }
        if (obj.nodeSubset.size() > 1) {
          IloRange cut;
          ctValue1 = 0.0;
          cut = obj.generatePREC(ctValue1, expr);
          expr.clear();
          if(ctValue1 > obj.nodeSubset.size() - 2 + EPSILON){
            add(cut, IloCplex::UseCutPurge).end();
            precFound = true;
            precCounter++;
          }
        }
      }
    }
    //cut 2: predecessor and successor cuts
    for (int i = 0; i < n; i++) {
      u = obj.pickups[i]; //pickup
      v = obj.sibling[u]; //correspondent delivery
      s = 2 * n + 2; //dummy node...
      t = u;
      Dinic graph = suppGraph;
      //modify the support graph: add the dummy node and the arcs going out of it...
      graph.addNode(); //V=2*n+3 0,1,...,2n+1,2n+2.
      for (int k = 0; k <= 2 * n + 1; k++) {
        int kap = 0;
        if (v != k && (obj.is_pickup[v] || obj.sibling[v] != k) && k != 0 && v != 2 * n + 1 && (!obj.is_pickup[v] || k != 2 * n + 1) && (v != 0 || obj.is_pickup[k])) {
          kap += (int) (10000 * obj._x[ obj.mapVarsX[100 * v + k] ] + 0.5);
        }
        if (u != k && (obj.is_pickup[k] || obj.sibling[k] != u) && u != 0 && k != 2 * n + 1 && (!obj.is_pickup[k] || u != 2 * n + 1) && (k != 0 || obj.is_pickup[u])) {
          kap += (int) (10000 * obj._x[ obj.mapVarsX[100 * k + u] ] + 0.5);
        }
        if (kap > 0) {
          graph.addArc(2 * n + 2, k, kap);
        }
      }
      graph.addArc(2 * n + 1, u, 20000); //this arc does not exist on the support graph
      flag = graph.modifyArc(0, u, 20000); //first check if 0->u already exists
      if (!flag) {
        graph.addArc(0, u, 20000);
      }
      graph.addArc(v, u, 20000); //this arc also does not exist on the support graph
      flow = graph.maxFlow(s, t);
      if (flow < 19990) {
        graph.runMinCut();
        obj.nodeSubset.clear();
        for (int ii = 0; ii <= 2 * n + 1; ii++) {
          if (graph.isInMinCut(ii)) {
            obj.nodeSubset.push_back(ii);
          }
        }
        if (obj.nodeSubset.size() >= 1) {
          IloRange cut2, cut3;
          ctValue1 = ctValue2 = ctValue3 = 0.0;
          expr.clear();
          obj.nodeSubset.push_back(v);
          cut2 = obj.generatePredecessorCut(ctValue2, expr);
          expr.clear();
          size = obj.nodeSubset.size();
          obj.nodeSubset[size - 1] = u;
          cut3 = obj.generateSuccessorCut(ctValue3, expr);
          expr.clear();
          if (ctValue2 > size - 1 + EPSILON) {
            add(cut2, IloCplex::UseCutPurge).end();
            precFound = true;
            precCounter++;
          }
          if (ctValue3 > size - 1 + EPSILON) {
            add(cut3, IloCplex::UseCutPurge).end();
            precFound = true;
            precCounter++;
          }
        }
      }
      graph.runMinCut();
      obj.nodeSubset.clear();
      for (int ii = 1; ii < 2 * n + 1; ii++) {
        if (graph.isInMinCut(ii)) {
          obj.nodeSubset.push_back(ii);
        }
      }

      double qS = 0;
      flag = false;
      for (int ii = 0; ii < obj.nodeSubset.size(); ii++) {
        qS += obj.itemLength[obj.nodeSubset[ii]];
        if (obj.nodeSubset[ii] == 0 || obj.nodeSubset[ii] == 2 * n + 1) {
          flag = true;
        }
      }
      if (qS < 0) {
        qS = -qS;
      }
      if (qS > obj.Q && !flag) {
        ctValue1 = 0.0;
        IloRange cut = obj.generateRoundedCapacity(ctValue1, expr);
        if (ctValue1 + EPSILON < ceil(qS / (obj.K * obj.Q))) {
          add(cut, IloCplex::UseCutPurge).end();
          rCapCounter++;
        }
        expr.clear();
      }
    }

  bool violatedLifos[2 * n + 2][2 * n + 2];
  memset(violatedLifos, false, sizeof (violatedLifos));
  int pi, pj, di, dj;

#ifdef SEPARATE_LIFO_EXACTLY_IN_ILOUSER

  for (int i = 0; i < obj.pickups.size(); i++) {
    for (int j = 0; j < obj.pickups.size(); j++) {
      if (obj.pickups[i] == obj.pickups[j])
        continue;
      pi = obj.pickups[i];
      pj = obj.pickups[j];
      di = obj.sibling[pi];
      dj = obj.sibling[pj];
      u = di; // n+i assume o papel de coleta
      v = pj; // j assume o papel de entrega
      s = 2 * n + 2; //dummy node...
      t = u;
      Dinic graph = suppGraph;
      Dinic auxGraph1 = graph;
      graph.addNode(); //V=2*n+3 0,1,...,2n+1,2n+2.
      for (int k = 0; k <= 2 * n + 1; k++) {
        int kap = 0;
        //if (v != k && !(!obj.is_pickup[v] && obj.sibling[v] == k) && k != 0 && v != 2 * n + 1 && !(obj.is_pickup[v] && k == 2 * n + 1) && !(v == 0 && !obj.is_pickup[k])) {
        //  kap += (int) (10000 * obj._x[ obj.mapVarsX[100 * v + k] ] + 0.5);
        //}
        if (u != k && !(!obj.is_pickup[k] && obj.sibling[k] == u) && u != 0 && k != 2 * n + 1 && !(obj.is_pickup[k] && u == 2 * n + 1) && !(k == 0 && !obj.is_pickup[u])) {
          kap += (int) (10000 * obj._x[ obj.mapVarsX[100 * k + u] ] + 0.5);
        }
        if (kap > 0) {
          graph.addArc(2 * n + 2, k, kap);
        }
      }

      graph.addNode(); //V=2*n+3 0,1,...,2n+1,2n+2.
      for (int ii = 1; ii <= 2 * n; ii++) {
        if (ii != pj && !(!obj.is_pickup[ii] && obj.sibling[ii] == pj)) {
          if (obj._x[ obj.mapVarsX[100 * ii + pj] ] > 0.0001) {
              graph.addArc(ii, 2 * n + 3, (int) (10000 * obj._x[ obj.mapVarsX[100 * ii + pj] ] + 0.5));
          }
        }
      }

      graph.addArc(2*n+2,pj,20000);

      graph.addArc(2 * n + 1, t, 20000);

      graph.addArc(0, t, 20000);

      flag = graph.modifyArc(2*n+3, t, 20000);
      if (!flag)
        graph.addArc(2*n+3, t, 20000);

      flag = graph.modifyArc(dj, t, 20000);
      if (!flag)
        graph.addArc(dj, t, 20000);

      flag = graph.modifyArc(pi, t, 20000);
      if (!flag)
        graph.addArc(pi, t, 20000);

      Dinic auxGraph2 = graph;
      flow = graph.maxFlow(s, t);

      if (flow < 19990) {
        graph.runMinCut();
        obj.nodeSubset.clear();
        for (int ii = 0; ii <= 2 * n + 1; ii++) {
          if (graph.isInMinCut(ii)) {
            obj.nodeSubset.push_back(ii);
          }
        }

        double sum_i_l;
        double sum_j_l;
        for (int k = 0; k < obj.K; k++) {
          sum_i_l = 0;
          sum_j_l = 0;
          for (int l = 0; l < obj.K; l++) {
            if (l != k) {
              sum_i_l += 10000 * obj._y[obj.mapVarsY[10 * pi + l]];
              sum_j_l += 10000 * obj._y[obj.mapVarsY[10 * pj + l]];
            }
          }
          if (flow + sum_i_l + sum_j_l < 19990) {
            ctValue1 = 0.0;
            IloRange cut;
            if(lifted_lifo){
              cut = obj.generate_lifted_lifo_loading_inside(obj.pickups[i], obj.pickups[j], k, ctValue1, expr);
            } else {
              cut = obj.generateLIFOS(obj.pickups[i], obj.pickups[j], k, ctValue1, expr);
            }
            expr.clear();
            if (ctValue1 > obj.nodeSubset.size() + 2 + EPSILON) {
              violatedLifos[i][j] = true;
              add(cut, IloCplex::UseCutPurge);
              lifoCounter++;
            }
          }
        }
      }
    }
  }
#endif

#ifdef FIND_VIOLATED_PATH_CAPACITY_AND_LIFO_IN_USERCALLBACK

    bool flag_2n1, deadend;
    bool visited[MAXNODES];
    bool inS[MAXNODES];
    int path[MAXNODES]; //defines the path seobj.quence
    int pIndex; //for indexing the discovered path
    int component[MAXNODES];
    int used_stack[MAXNODES];
    vector< vector<int> > paths;
    int comp;
    memset(component, 0, sizeof (component));
    memset(path, 0, sizeof (path));
    memset(visited, false, sizeof (visited));

    /* Used stacks */
    for (int i = 0; i < obj.pickups.size(); i++) {
      for (int k = 0; k < obj.K; k++) {
        if (obj._y[ obj.mapVarsY[10 * obj.pickups[i] + k] ] > 0.1) {
          used_stack[obj.pickups[i] ] = k;
          used_stack[obj.sibling[obj.pickups[i] ]] = k;
          k = obj.K;
        }
      }
    }

    /* Path starting at 0 */
    u = 0;
    pIndex = 0;
    comp = 0;
    vector<int> current_path;
    current_path.push_back(u);
    flag_2n1 = false;
    deadend = false;
    while (!flag_2n1 && !deadend) {
      deadend = true;
      for (int v = 1; v <= 2 * n + 1; v++) {
        if (u != v && !(!obj.is_pickup[u] && obj.sibling[u] == v) && u != 2 * n + 1
                && !(obj.is_pickup[u] && v == 2 * n + 1) && !(u == 0 && !obj.is_pickup[v])
                && !(u == 0 && v == (2 * n + 1))) {
          if (obj._x[ obj.mapVarsX[100 * u + v] ] > 0.8) {
            deadend = false;
            component[v] = comp;
            visited[v] = true;
            u = v;
            v = 2 * n + 2;
            if (u == 2 * n + 1) {
              flag_2n1 = true;
            }
            current_path.push_back(u);
          }
        }
      }
    }
    if (current_path.size() > obj.Q) {
      paths.push_back(current_path);
    }

    for (int i = 0; i < obj.pickups.size(); i++) {
      if (!visited[obj.pickups[i]]) {
        /* Forward Path */
        current_path.clear();
        flag_2n1 = false;
        deadend = false;
        comp++;
        u = obj.pickups[i];
        while (!flag_2n1 && !deadend) {
          deadend = true;
          for (int v = 1; v <= 2 * n + 1; v++) {
            if (u != v && !(!obj.is_pickup[u] && obj.sibling[u] == v) && u != 2 * n + 1
                    && !(obj.is_pickup[u] && v == 2 * n + 1) && !(u == 0 && !obj.is_pickup[v])) {
              if (obj._x[ obj.mapVarsX[100 * u + v] ] > 0.8) {
                if (visited[v]) { /* If v was already visited it is cycle */
                  deadend = true;
                } else {
                  deadend = false;
                  component[v] = comp;
                  visited[v] = true;
                  u = v;
                  v = 2 * n + 2;
                  if (u == 2 * n + 1)
                    flag_2n1 = true;
                  current_path.push_back(u);
                }
              }
            }
          }
        }
        if (current_path.size() > obj.Q) {
          paths.push_back(current_path);
        }
      }
    }

    //lifoFound = false;
    //conflictFound = false;
    for (int p = 0; p < paths.size(); p++) {


      short int start[MAXSTKS];
      short int end[MAXSTKS];
      bool visited2[MAXNODES];
      short int capacities[MAXSTKS];
      short int old_capacity;
      vector<int> loaded_in_k;
      obj.nodeSubset.clear();
      pIndex = paths[p].size() - 1;

      u = paths[p][0];
      for (int i = 1; i < pIndex; i++) {
        v = paths[p][i];
        if (obj.is_pickup[v] && component[obj.sibling[v]] == 1) //v is a pickup node
        {
          short int stk = used_stack[v];
          short int dv = obj.sibling[v];
          short int last = 0;
          memset(visited, false, sizeof (visited));
          u = v;
          for (int ii = i + 1; ii < pIndex; ii++) {
            v = paths[p][ii];
            if (v != dv) {
              if (used_stack[v] == stk) {
                visited[v] = true;
              }
              u = v;
            } else {
              last = ii;
              ii = pIndex;
            }
          }

          flag = false;
          obj.nodeSubset.clear();

          if (lifted_lifo == true) {

            bool possible_lifo = false;
            vector<int> loaded_in_k;
            memset(inS, false, sizeof (inS));

            for (int j = last - 1; j > i; --j) {
              u = paths[p][j];
              v = obj.sibling[u];
              obj.nodeSubset.push_back(u);
              inS[u] = true;
              if (obj.is_pickup[u] && visited[u] && !visited[v]) {
                loaded_in_k.push_back(u);
                j = i;
                possible_lifo = true;
              }
            }

            if (possible_lifo && !violatedLifos[paths[p][i]][loaded_in_k[0]]) {
              for (int j = 0; j < obj.nodeSubset.size(); j++) {
                u = obj.nodeSubset[j];
                v = obj.sibling[u];
                if (obj.is_pickup[u] && !inS[v] && obj.nodeSubset[j] != loaded_in_k[0]) {
                  loaded_in_k.push_back(u);
                }
              }
              for (int kk = 0; kk < obj.K; kk++) {
                ctValue1 = 0.0;
                IloRange cut = obj.generate_lifted_lifo(paths[p][i], loaded_in_k[0], kk, ctValue1, expr, loaded_in_k);
                expr.clear();
                if (ctValue1 > obj.nodeSubset.size() + 2 + EPSILON) {
                  violatedLifos[paths[p][i]][loaded_in_k[0]] = true;
                  add(cut, IloCplex::UseCutPurge);
                  lifoCounter++;
                  flag = true;
                }
              }
            }

            loaded_in_k.clear();
            obj.nodeSubset.clear();
            possible_lifo = false;
            memset(inS, false, sizeof (inS));

            obj.nodeSubset.push_back(paths[p][i]);
            loaded_in_k.push_back(paths[p][i]);

            int last_delivery = -1;
            for (int j = i + 1; j < last; ++j) {
              u = paths[p][j];
              v = obj.sibling[u];
              obj.nodeSubset.push_back(u);
              inS[u] = true;
              if (!obj.is_pickup[u] && visited[u] && !visited[v]) {
                last_delivery = j;
                j = last;
                possible_lifo = true;
                obj.nodeSubset.pop_back();
              }
            }

            if (possible_lifo && violatedLifos[obj.sibling[paths[p][last_delivery]]][paths[p][i]]) {
              for (int j = 0; j < obj.nodeSubset.size(); j++) {
                u = obj.nodeSubset[j];
                v = obj.sibling[u];
                if (obj.is_pickup[u] && !inS[v] && obj.nodeSubset[j] != loaded_in_k[0]) {
                  loaded_in_k.push_back(u);
                }
              }
              for (int kk = 0; kk < obj.K; kk++) {
                ctValue1 = 0.0;
                IloRange cut = obj.generate_lifted_lifo(obj.sibling[paths[p][last_delivery]], paths[p][i], kk, ctValue1, expr, loaded_in_k);
                expr.clear();
                if (ctValue1 > obj.nodeSubset.size() + 2 + EPSILON) {
                  violatedLifos[obj.sibling[paths[p][last_delivery]]][paths[p][i]] = true;
                  add(cut, IloCplex::UseCutPurge);
                  lifoCounter++;
                  flag = true;
                }
              }
            }

            loaded_in_k.clear();
            obj.nodeSubset.clear();

          }


          if (original_lifo == true) {
            for (int j = last - 1; j >= i + 1; j--) {
              u = paths[p][j];
              v = obj.sibling[u];
              obj.nodeSubset.push_back(u);
              if (obj.is_pickup[u]) {
                if (visited[u] && !visited[v])//a pickup in S for which the delivery is outside S
                {
                  for (int kk = 0; kk < obj.K; kk++) {
                    ctValue1 = 0.0;
                    IloRange cut = obj.generateLIFOS(paths[p][i], u, kk, ctValue1, expr);
                    expr.clear();
                    if (ctValue1 > obj.nodeSubset.size() + 2 + EPSILON) {
                      add(cut, IloCplex::UseCutPurge);
                      lifoCounter++;
                      flag = true;
                      j = 0;
                    }
                  }
                }
              }
            }

            obj.nodeSubset.clear();

            obj.nodeSubset.push_back(paths[p][i]);

            for (int j = i + 1; j < last; j++) {
              u = paths[p][j];
              v = obj.sibling[u];
              obj.nodeSubset.push_back(u);
              if (!obj.is_pickup[u]) {
                if (visited[u] && !visited[v])//a delivery in S for which the pickup is outside S
                {
                  visited[u] = false;
                  obj.nodeSubset.pop_back();
                  for (int kk = 0; kk < obj.K; kk++) {
                    ctValue1 = 0.0;
                    IloRange cut = obj.generateLIFOS(v, paths[p][i], kk, ctValue1, expr);
                    expr.clear();
                    if (ctValue1 > obj.nodeSubset.size() + 2 + EPSILON) {
                      add(cut, IloCplex::UseCutPurge);
                      lifoCounter++;
                      flag = true;
                      j = last;
                    }
                  }
                }
              }
            }

            obj.nodeSubset.clear();

          }
        }

      }

      expr.clear();
      memset(capacities, 0, sizeof (capacities));
      memset(visited2, false, sizeof (visited2));
      
      obj.nodeSubset.clear();
      
      if (sec_capacity) {
        for (int i = 1; i < pIndex; i++) {
          v = paths[p][i];
          int k = used_stack[v];
          if (obj.is_pickup[v]) {
            old_capacity = capacities[k];
            capacities[k] += obj.itemLength[v];
            visited2[v] = true;
            if (old_capacity == 0) {//before loading v, stack k was empty...
              start[k] = i - 1;
            } else if (capacities[k] > obj.Q) {
              end[k] = i;
              obj.nodeSubset.push_back(paths[p][start[k]]);
              for (int ti = start[k]; ti < end[k]; ti++) {
                if (obj.is_pickup[paths[p][ti + 1]] && used_stack[paths[p][ti + 1]] == k) {
                  obj.nodeSubset.push_back(paths[p][ti + 1]);
                  loaded_in_k.push_back(paths[p][ti + 1]);
                } else {
                  obj.nodeSubset.push_back(paths[p][ti + 1]);
                }
              }

              for (int kk = 0; kk < obj.K; kk++) {
                ctValue1 = 0.0;
                IloRange cut = obj.generateSECCapacity(kk, ctValue1, expr, loaded_in_k);
                expr.clear();
                if (ctValue1 > loaded_in_k.size() + obj.nodeSubset.size() - 2 + EPSILON) {
                  add(cut, IloCplex::UseCutPurge);
                  pCapCounter++;
                }
              }

              loaded_in_k.clear();
              obj.nodeSubset.clear();
            }
          }
          if (!obj.is_pickup[v] && visited2[obj.sibling[v]]) {
            capacities[k] -= -obj.itemLength[v];
          }
          u = v;
        }
      }

      if (path_capacity) {

        int count = 0;
        ctValue1 = 0;
        for (int i = 1; i < pIndex; i++) {
          v = paths[p][i];
          int k = used_stack[v];
          if (obj.is_pickup[v]) {
            old_capacity = capacities[k];
            capacities[k] += obj.itemLength[v];
            visited2[v] = true;
            if (old_capacity == 0) //before loading v, stack k was empty...
            {
              start[k] = i - 1;
            } else if (capacities[k] > obj.Q) {
              end[k] = i;
              for (int kk = 0; kk < obj.K; kk++) {//generate a paths[p] inequality
                for (int ti = start[k]; ti < end[k]; ti++) {
                  if (obj.is_pickup[paths[p][ti + 1]] && used_stack[paths[p][ti + 1]] == k) {
                    expr += obj.y[ obj.mapVarsY[10 * paths[p][ti + 1] + kk] ];
                    expr += obj.x[ obj.mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ];
                    ctValue1 += obj._x[ obj.mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ] + obj._y[ obj.mapVarsY[10 * paths[p][ti + 1] + kk] ];
                    count += 2;
                  } else {
                    count++;
                    expr += obj.x[ obj.mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ];
                    ctValue1 += obj._x[ obj.mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ];
                  }
                }
                if (ctValue1 > count - 1 + EPSILON) {
                  add(expr <= count - 1, IloCplex::UseCutPurge);
                  pCapCounter++;
                }
                expr.clear();
                count = 0;
                ctValue1 = 0;
              }
            }
          }
          if (!obj.is_pickup[v] && visited2[obj.sibling[v]]) {
            capacities[k] -= -obj.itemLength[v];
          }
        }
      }

    }

#endif

  } catch (IloException& ex) {
    cout << "userCallback: " << ex << endl;
  }
}

ILOINCUMBENTCALLBACK1(incumbentCallback, TSPPDMS&, obj) {
  obj.setUB(getObjValue());
}

void TSPPDMS::setUB(IloNum value) {
  ub = value;
}

IloRange TSPPDMS::search_violated_rounded_capacity_inequalities(bool &foundCut) {
  nodeSubset.clear();
  int starting_node;
  double lambda1, lambda2, lambda3;
  double exp1, exp2, exp3;
  double lb_entering, lb_leaving;
  double xDeltaPlusS;
  bool inS[2 * n + 1];
  memset(inS, false, sizeof (inS));
  lb_leaving = lb_entering = 0;
  lambda1 = 5 * ((double) rand() / (RAND_MAX));
  lambda2 = 5 * ((double) rand() / (RAND_MAX));
  lambda3 = (double) rand() / (RAND_MAX);
  starting_node = rand() % (2 * n) + 1; // [1 ... 2n]
  nodeSubset.push_back(starting_node);
  inS[starting_node] = true;
  if (is_pickup[starting_node]) {
    lb_entering = 0;
    lb_leaving += itemLength[sibling[starting_node]];
  } else {
    lb_entering += itemLength[sibling[starting_node]];
    lb_leaving = 0;
  }
  xDeltaPlusS = 1.0;
  bool stop_condition = false;
  double greatest_increase = 0.0;
  int greatest_increase_node;
  while (!stop_condition) {
    greatest_increase = -50000;
    greatest_increase_node = -1;
    for (int i = 1; i < 2 * n + 1; i++) {
      if (!inS[i]) {
        if (is_pickup[i]) {
          if (inS[sibling[i]]) {
            lb_entering -= itemLength[i];
          } else {
            lb_leaving += itemLength[sibling[i]];
          }
        } else {
          if (inS[sibling[i]]) {
            lb_leaving -= itemLength[i];
          } else {
            lb_entering += itemLength[sibling[i]];
          }
        }
        nodeSubset.push_back(i); // Add to S to make the following sum easier
        inS[i] = true;
        xDeltaPlusS = 0.0;
        for (int s = 0; s < nodeSubset.size(); s++) {
          int u = nodeSubset[s];
          for (int j = 1; j <= 2 * n + 1; j++) {
            if (!inS[j] && u != j && !(!is_pickup[u] && sibling[u] == j) && !(is_pickup[u] && j == (2 * n + 1))) {
              xDeltaPlusS += _x[mapVarsX[100 * nodeSubset[s] + j]];
            }
          }
        }
        exp1 = lambda1 * (max(lb_entering, -1 * lb_leaving) - Q * xDeltaPlusS);
        exp2 = lambda2 * Q * (max(ceil(lb_entering / Q), ceil(-1 * lb_leaving / Q)) - xDeltaPlusS);
        exp3 = lambda3 * (min(lb_entering, -1 * lb_leaving) - Q * xDeltaPlusS);
        if (exp1 + exp2 + exp3 > greatest_increase) {
          greatest_increase = exp1 + exp2 + exp3;
          greatest_increase_node = i;
        }

        nodeSubset.pop_back();
        inS[i] = false;
        if (is_pickup[i]) {
          if (inS[sibling[i]]) {
            lb_entering += itemLength[i];
          } else {
            lb_leaving -= itemLength[sibling[i]];
          }
        } else {
          if (inS[sibling[i]]) {
            lb_leaving += itemLength[i];
          } else {
            lb_entering -= itemLength[sibling[i]];
          }
        }

      }
    }
    if (greatest_increase_node != -1) {
      //      printf("Added node %d to S\n", greatest_increase_node);
      nodeSubset.push_back(greatest_increase_node);
      inS[greatest_increase_node] = true;
      if (is_pickup[greatest_increase_node]) {
        if (inS[sibling[greatest_increase_node]]) {
          lb_entering -= itemLength[greatest_increase_node];
        } else {
          lb_leaving += itemLength[sibling[greatest_increase_node]];
        }
      } else {
        if (inS[sibling[greatest_increase_node]]) {
          lb_leaving -= itemLength[greatest_increase_node];
        } else {
          lb_entering += itemLength[sibling[greatest_increase_node]];
        }
      }

      if ((max(ceil(lb_entering / Q), ceil(-1 * lb_leaving) / Q) - xDeltaPlusS) > EPSILON) {
        stop_condition = true;
      }
    } else { // Didn't happen any increase in this iteration
      stop_condition = true;
    }
  }

  double ctValue1 = 0.0;
  IloExpr expr(Env);
  IloRange cut = generateRoundedCapacity(ctValue1, expr);
  double qS = max(lb_entering,-1*lb_leaving);

  if (ctValue1 + EPSILON < ceil(qS / (K * Q))) {
    foundCut = true;
  }
  return cut;
  nodeSubset.clear();
}

bool TSPPDMS::findCuts(int& counter1, int& counter2, int& counter3, int& counter4) {
  int s, t, u, v;
  int flow, size;
  bool incS[MAXNODES];
  IloExpr expr(Env);
  IloNum ctValue1, ctValue2, ctValue3;
  bool flag;
  bool precFound, capFound;
  fcuts++;
  //printf("\n### Findcuts %d ###\n", fcuts);
  //  create_graphviz_image(fcuts, LINEAR_SOLUTION);

  for (int i = 0; i < 5; i++) {
    bool foundCut = false;
    IloRange cut = search_violated_rounded_capacity_inequalities(foundCut);
    if (foundCut) {
      caps.add(cut);
      rCapCounter++;
      i=5;
    }
  }

  Dinic suppGraph(2 * n + 2);
  for (int i = 1; i <= 2 * n; i++) {
    for (int j = 1; j <= 2 * n; j++) {
      if (i != j && !(!is_pickup[i] && sibling[i] == j)) {
        if (_x[ mapVarsX[100 * i + j] ] > 0.0001) {
          suppGraph.addArc(i, j, (int) (10000 * _x[ mapVarsX[100 * i + j] ] + 0.5));
        }
      }
    }
  }

  for (int i = 0; i < pickups.size(); i++) {
    if (_x[ mapVarsX[100 * 0 + pickups[i]] ] > 0.0001) {
      suppGraph.addArc(0, pickups[i], (int) (10000 * _x[ mapVarsX[100 * 0 + pickups[i]] ] + 0.5));
    }
  }

  for (int i = 0; i < deliveries.size(); i++) {
    if (_x[ mapVarsX[100 * deliveries[i]+(2 * n + 1)] ] > 0.0001) {
      suppGraph.addArc(deliveries[i], 2 * n + 1, (int) (10000 * _x[ mapVarsX[100 * deliveries[i]+(2 * n + 1)] ] + 0.5));
    } 
  }

  secs.clear();
  precs.clear();
  lifos.clear();
  caps.clear();
  /****************************************************************/
  /*************SUBTOUR ELIMINATION CONSTRAINTS********************/
  /****************************************************************/
  //solves the maxflow from depot 0 to each of the pickups
  s = 0;
  memset(incS, false, sizeof (incS));
  capFound = false;
  for (int i = 0; i < pickups.size(); i++) {
    t = pickups[i];
    if (!incS[t]) {
      Dinic graph = suppGraph;
      flow = graph.maxFlow(s, t);
      if (flow < 9990) {

        graph.runMinCut();
        nodeSubset.clear();
        for (int ii = 0; ii <= 2 * n + 1; ii++) {
          if (!graph.isInMinCut(ii)) {
            nodeSubset.push_back(ii);
          }
        }

        if (nodeSubset.size() > 1) {
          IloRange cut;
          ctValue1 = 0.0;
          cut = generateSEC(ctValue1, expr);
          expr.clear();
          if(ctValue1 > nodeSubset.size() - 1 + EPSILON){
            secs.add(cut);
            secCounter++;            
            for (int ii = 0; ii < nodeSubset.size(); ii++) {
              incS[nodeSubset[ii]] = true;
            }
          }
        }
      } else if (!capFound) {
        graph.runMinCut();
        nodeSubset.clear();
        for (int ii = 1; ii < 2 * n + 1; ii++) {
          if (graph.isInMinCut(ii)) {
            nodeSubset.push_back(ii);
          }
        }
        double qS = 0;
        flag = false;
        for (int ii = 0; ii < nodeSubset.size(); ii++) {
          qS += itemLength[nodeSubset[ii]];
          if (nodeSubset[ii] == 0 || nodeSubset[ii] == 2 * n + 1) {
            flag = true;
          }
        }
        if (qS < 0) {
          qS = -qS;
        }
        if (qS > Q && !flag) {
          ctValue1 = 0.0;
          IloRange cut = generateRoundedCapacity(ctValue1, expr);
          if (ctValue1 + EPSILON < ceil(qS / (K * Q))) {
            caps.add(cut);
            capFound = true;
            rCapCounter++;
          }
          expr.clear();
        }
      }
    }
  }
  //solves the maxflow from each of the deliveires to the depot 2n+1  t = 2 * n + 1;
  memset(incS, false, sizeof (incS));
  capFound = false;
  //printf("To 2n+1\n");
  for (int i = 0; i < deliveries.size(); i++) {
    s = deliveries[i];
    if (!incS[s]) {
      Dinic graph = suppGraph;
      flow = graph.maxFlow(s, t);
      if (flow < 9990) {
        graph.runMinCut();
        nodeSubset.clear();
        for (int ii = 0; ii <= 2 * n + 1; ii++) {
          if (graph.isInMinCut(ii)) {
            nodeSubset.push_back(ii);
          }
        }
        if (nodeSubset.size() > 1) {
          IloRange cut;
          ctValue1 = 0.0;
          cut = generateSEC(ctValue1, expr);
          expr.clear();
          if(ctValue1 > nodeSubset.size() - 1 + EPSILON){
            secs.add(cut);
            secCounter++;            
            for (int ii = 0; ii < nodeSubset.size(); ii++) {
              incS[nodeSubset[ii]] = true;
            }
          }
        }
      } else if (!capFound) {
        graph.runMinCut();
        nodeSubset.clear();
        for (int ii = 1; ii < 2 * n + 1; ii++) {
          if (graph.isInMinCut(ii)) {
            nodeSubset.push_back(ii);
          }
        }
        double qS = 0;
        flag = false;
        for (int ii = 0; ii < nodeSubset.size(); ii++) {
          qS += itemLength[nodeSubset[ii]];
          if (nodeSubset[ii] == 0 || nodeSubset[ii] == 2 * n + 1) {
            flag = true;
          }
        }
        if (qS < 0) {
          qS = -qS;
        }
        if (qS > Q && !flag) {
          ctValue1 = 0.0;
          IloRange cut = generateRoundedCapacity(ctValue1, expr);
          if (ctValue1 + EPSILON < ceil(qS / (K * Q))) {
            caps.add(cut);
            capFound = true;
            rCapCounter++;
          }
          expr.clear();
        }
      }
    }
  }

  /****************************************************************/
  /*************PRECEDENCE CONSTRAINTS*****************************/
  /****************************************************************/
  precFound = false;
  for (int i = 0; i < n; i++) {
    u = pickups[i];
    v = sibling[u];
    s = 0;
    t = 2 * n + 1;
    Dinic graph = suppGraph;
    graph.addArc(s, v, 20000);
    graph.addArc(u, t, 20000);
    flow = graph.maxFlow(s, t);
    if (flow < 19990) {
      graph.runMinCut();
      nodeSubset.clear();
      for (int ii = 0; ii <= 2 * n + 1; ii++) {
        if (graph.isInMinCut(ii)) {
          nodeSubset.push_back(ii);
        }
      }
      if (nodeSubset.size() > 1) {
        IloRange cut;
        ctValue1 = 0.0;
        cut = generatePREC(ctValue1, expr);
        expr.clear();
        if(ctValue1 > nodeSubset.size() - 2 + EPSILON){
          precs.add(cut);
          precFound = true;
          precCounter++;
        }
      }
    }
  }
  //cut 2: predecessor and successor cuts
  for (int i = 0; i < n; i++) {
    u = pickups[i]; //pickup
    v = sibling[u]; //correspondent delivery
    s = 2 * n + 2; //dummy node...
    t = u;
    Dinic graph = suppGraph;
    //modify the support graph: add the dummy node and the arcs going out of it...
    graph.addNode(); //V=2*n+3 0,1,...,2n+1,2n+2.
    for (int k = 0; k <= 2 * n + 1; k++) {
      //cap[2*d.n+2][k] = cap[v][k] + cap[k][u];
      int kap = 0;
      if (v != k && (is_pickup[v] || sibling[v] != k) && k != 0 && v != 2 * n + 1 && (!is_pickup[v] || k != 2 * n + 1) && (v != 0 || is_pickup[k])) {
        kap += (int) (10000 * _x[ mapVarsX[100 * v + k] ] + 0.5);
      }
      if (u != k && (is_pickup[k] || sibling[k] != u) && u != 0 && k != 2 * n + 1 && (!is_pickup[k] || u != 2 * n + 1) && (k != 0 || is_pickup[u])) {
        kap += (int) (10000 * _x[ mapVarsX[100 * k + u] ] + 0.5);
      }
      if (kap > 0) {
        graph.addArc(2 * n + 2, k, kap);
      }
    }
    graph.addArc(2 * n + 1, t, 20000); //this arc does not exist on the support graph
    flag = graph.modifyArc(0, t, 20000); //first check if 0->u already exists
    if (!flag)
      graph.addArc(0, t, 20000);
    graph.addArc(v, t, 20000); //this arc also does not exist on the support graph
    flow = graph.maxFlow(s, t);
    if (flow < 19990) {
      graph.runMinCut();
      nodeSubset.clear();
      for (int ii = 0; ii <= 2 * n + 1; ii++)
        if (graph.isInMinCut(ii)) {
          nodeSubset.push_back(ii);
        }
      if (nodeSubset.size() >= 2) {
        IloRange cut2, cut3;
        ctValue1 = ctValue2 = ctValue3 = 0.0;
        //        cut1 = generatePRECI(ctValue1, expr, u); //x(n+i,S)+x(S)+x(S,i)<=|S|
        expr.clear();
        nodeSubset.push_back(v);
        cut2 = generatePredecessorCut(ctValue2, expr);
        expr.clear();
        size = nodeSubset.size();
        nodeSubset[size - 1] = u;
        cut3 = generateSuccessorCut(ctValue3, expr);
        expr.clear();
        if (ctValue2 > size - 1)//violated predecessor inequality
        {
          precs.add(cut2);
          npredCuts++;
          precCounter++;
        }
        if (ctValue3 > size - 1)//violated successor inequality
        {
          precs.add(cut3);
          nsuccCuts++;
          precCounter++;
        }
      }
    }

    graph.runMinCut();
    nodeSubset.clear();
    for (int ii = 1; ii < 2 * n + 1; ii++) {
      if (graph.isInMinCut(ii)) {
        nodeSubset.push_back(ii);
      }
    }

    double qS = 0;
    flag = false;
    for (int ii = 0; ii < nodeSubset.size(); ii++) {
      qS += itemLength[nodeSubset[ii]];
      if (nodeSubset[ii] == 0 || nodeSubset[ii] == 2 * n + 1) {
        flag = true;
      }
    }
    if (qS < 0) {
      qS = -qS;
    }
    if (qS > Q && !flag) {
      ctValue1 = 0.0;
      IloRange cut = generateRoundedCapacity(ctValue1, expr);
      if (ctValue1 + EPSILON < ceil(qS / (K * Q))) {
        caps.add(cut);
        rCapCounter++;
      }
      expr.clear();
    }
  }

  nodeSubset.clear();

  bool violatedLifos[2 * n + 2][2 * n + 2];
  memset(violatedLifos, false, sizeof (violatedLifos));
  int pi, pj, di, dj;

  for (int i = 0; i < pickups.size(); i++) {
    for (int j = 0; j < pickups.size(); j++) {
      if (pickups[i] == pickups[j])
        continue;
      pi = pickups[i];
      pj = pickups[j];
      di = sibling[pi];
      dj = sibling[pj];
      u = di; // n+i assume o papel de coleta
      v = pj; // j assume o papel de entrega
      s = 2 * n + 2; //dummy node...
      t = u;
      Dinic graph = suppGraph;
      Dinic auxGraph1 = graph;
      //modify the support graph: add the dummy node and the arcs going out of it...
      graph.addNode(); //V=2*n+3 0,1,...,2n+1,2n+2.
      for (int k = 0; k <= 2 * n + 1; k++) {
        int kap = 0;
        //if (v != k && !(!is_pickup[v] && sibling[v] == k) && k != 0 && v != 2 * n + 1 && !(is_pickup[v] && k == 2 * n + 1) && !(v == 0 && !is_pickup[k])) {
        //  kap += (int) (10000 * _x[ mapVarsX[100 * v + k] ] + 0.5);
        //}
        if (u != k && !(!is_pickup[k] && sibling[k] == u) && u != 0 && k != 2 * n + 1 && !(is_pickup[k] && u == 2 * n + 1) && !(k == 0 && !is_pickup[u])) {
          kap += (int) (10000 * _x[ mapVarsX[100 * k + u] ] + 0.5);
        }
        if (kap > 0) {
          graph.addArc(2 * n + 2, k, kap);
        }
      }

      graph.addNode(); //V=2*n+3 0,1,...,2n+1,2n+2.
      for (int ii = 1; ii <= 2 * n; ii++) {
        if (ii != pj && !(!is_pickup[ii] && sibling[ii] == pj)) {
          if (_x[ mapVarsX[100 * ii + pj] ] > 0.0001) {
              graph.addArc(ii, 2 * n + 3, (int) (10000 * _x[ mapVarsX[100 * ii + pj] ] + 0.5));
          }
        }
      }

      /*for (int k = 0; k <= 2 * n + 1; k++) {
        int kap = 0;
        if (v != k && !(!is_pickup[k] && sibling[k] == v) && v != 0 && k != 2 * n + 1 && !(is_pickup[k] && v == 2 * n + 1) && !(k == 0 && !is_pickup[v])) {
          kap += (int) (10000 * _x[ mapVarsX[100 * k + v] ] + 0.5);
        }
        if (kap > 0) {
          graph.addArc(2 * n + 2, k, kap);
        }
      }*/

      graph.addArc(2*n+2,pj,20000);

      graph.addArc(2 * n + 1, t, 20000);

      graph.addArc(0, t, 20000);

      flag = graph.modifyArc(2*n+3, t, 20000);
      if (!flag)
        graph.addArc(2*n+3, t, 20000);

      flag = graph.modifyArc(dj, t, 20000);
      if (!flag)
        graph.addArc(dj, t, 20000);

      flag = graph.modifyArc(pi, t, 20000);
      if (!flag)
        graph.addArc(pi, t, 20000);

      Dinic auxGraph2 = graph;
      flow = graph.maxFlow(s, t);

      if (flow < 19990) {
        graph.runMinCut();
        nodeSubset.clear();
        for (int ii = 0; ii <= 2 * n + 1; ii++) {
          if (graph.isInMinCut(ii)) {
            nodeSubset.push_back(ii);
          }
        }

        double sum_i_l;
        double sum_j_l;
        for (int k = 0; k < K; k++) {
          sum_i_l = 0;
          sum_j_l = 0;
          for (int l = 0; l < K; l++) {
            if (l != k) {
              sum_i_l += 10000 * _y[mapVarsY[10 * pi + l]];
              sum_j_l += 10000 * _y[mapVarsY[10 * pj + l]];
            }
          }
          if (flow + sum_i_l + sum_j_l < 19990) {
            ctValue1 = 0.0;
            IloRange cut;
            if(lifted_lifo){
              cut = generate_lifted_lifo_loading_inside(pickups[i], pickups[j], k, ctValue1, expr);
            } else {
              cut = generateLIFOS(pickups[i], pickups[j], k, ctValue1, expr);
            }
            expr.clear();
            if (ctValue1 > nodeSubset.size() + 2 + EPSILON) {
              violatedLifos[i][j] = true;
              lifos.add(cut);
              lifoCounter++;
              /*printf("Maxflow from %d to %d equals to %d",s,t,flow);
              printf("%d %d : S: ", pi, pj);
              for (int a = 0; a < nodeSubset.size(); a++) {
                printf("%d ", nodeSubset[a]);
              }
              printf("\n");
              create_graphviz_image(fcuts, LINEAR_SOLUTION);
              char filename[60];
              sprintf(filename, "%d_lifo_%d_between_%d_%d_support_graph", fcuts, lifoCounter, pickups[i], pickups[j]);
              auxGraph1.create_graphviz_image(filename, pickups, deliveries);
              sprintf(filename, "%d_lifo_%d_between_%d_%d_support_graph_modified", fcuts, lifoCounter, pickups[i], pickups[j]);
              auxGraph2.create_graphviz_image(filename, pickups, deliveries);*/
            }
          }
        }
      }
    }
  }

#ifdef FIND_VIOLATED_PATH_CAPACITY_AND_LIFO_IN_FINDCUTS

  bool flag_2n1, deadend;
  bool visited[MAXNODES];
  int path[MAXNODES]; //defines the path sequence
  bool inS[MAXNODES];
  int pIndex; //for indexing the discovered path
  int component[MAXNODES];
  int used_stack[MAXNODES];
  vector< vector<int> > paths;
  int comp;
  memset(component, 0, sizeof (component));
  memset(path, 0, sizeof (path));
  memset(visited, false, sizeof (visited));

  /* Used stacks */
  for (int i = 0; i < pickups.size(); i++) {
    for (int k = 0; k < K; k++) {
      if (_y[ mapVarsY[10 * pickups[i] + k] ] > 0.1) {
        used_stack[pickups[i] ] = k;
        used_stack[sibling[pickups[i] ]] = k;
        k = K;
      }
    }
  }

  /* Path starting at 0 */
  u = 0;
  pIndex = 0;
  comp = 0;
  vector<int> current_path;
  current_path.push_back(u);
  flag_2n1 = false;
  deadend = false;
  while (!flag_2n1 && !deadend) {
    deadend = true;
    for (int v = 1; v <= 2 * n + 1; v++) {
      if (u != v && !(!is_pickup[u] && sibling[u] == v) && u != 2 * n + 1
              && !(is_pickup[u] && v == 2 * n + 1) && !(u == 0 && !is_pickup[v])
              && !(u == 0 && v == (2 * n + 1))) {
        if (_x[ mapVarsX[100 * u + v] ] > 0.8) {
          deadend = false;
          component[v] = comp;
          visited[v] = true;
          u = v;
          v = 2 * n + 2;
          if (u == 2 * n + 1) {
            flag_2n1 = true;
          }
          current_path.push_back(u);
        }
      }
    }
  }
  if (current_path.size() > Q) {
    paths.push_back(current_path);
  }

  for (int i = 0; i < pickups.size(); i++) {
    if (!visited[pickups[i]]) {
      /* Forward Path */
      current_path.clear();
      flag_2n1 = false;
      deadend = false;
      comp++;
      u = pickups[i];
      while (!flag_2n1 && !deadend) {
        deadend = true;
        for (int v = 1; v <= 2 * n + 1; v++) {
          if (u != v && !(!is_pickup[u] && sibling[u] == v) && u != 2 * n + 1
                  && !(is_pickup[u] && v == 2 * n + 1) && !(u == 0 && !is_pickup[v])) {
            if (_x[ mapVarsX[100 * u + v] ] > 0.8) {
              if (visited[v]) { /* If v was already visited it is cycle */
                deadend = true;
              } else {
                deadend = false;
                component[v] = comp;
                visited[v] = true;
                u = v;
                v = 2 * n + 2;
                if (u == 2 * n + 1)
                  flag_2n1 = true;
                current_path.push_back(u);
              }
            }
          }
        }
      }
      if (current_path.size() > Q) {
        paths.push_back(current_path);
      }
    }
  }

  for (int p = 0; p < paths.size(); p++) {

    short int start[MAXSTKS];
    short int end[MAXSTKS];
    bool visited2[MAXNODES];
    short int capacities[MAXSTKS];
    short int old_capacity;
    vector<int> loaded_in_k;
    nodeSubset.clear();

    pIndex = paths[p].size() - 1;

    u = paths[p][0];
    for (int i = 1; i < pIndex; i++) {
      v = paths[p][i];
      if (is_pickup[v] && component[sibling[v]] == component[v]) //v is a pickup node
      {
        short int stk = used_stack[v];
        short int dv = sibling[v];
        short int last = 0;
        memset(visited, false, sizeof (visited));
        //the stack being checked

        u = v;
        for (int ii = i + 1; ii < pIndex; ii++) {
          v = paths[p][ii];
          if (v != dv) {
            if (used_stack[v] == stk) {
              visited[v] = true;
            }
            u = v;
          } else {
            last = ii;
            ii = pIndex; //breaks the for...
          }
        }

        flag = false;
        nodeSubset.clear();

        if (lifted_lifo) {

          memset(inS, false, sizeof (inS));
          bool possible_lifo = false;
          vector<int> loaded_in_k;

          for (int j = last - 1; j > i; --j) {
            u = paths[p][j];
            v = sibling[u];
            nodeSubset.push_back(u);
            inS[u] = true;
            if (is_pickup[u] && visited[u] && !visited[v]) {
              loaded_in_k.push_back(u);
              j = i;
              possible_lifo = true;
            }
          }

          if (possible_lifo && !violatedLifos[paths[p][i]][loaded_in_k[0]]) {
            for (int j = 0; j < nodeSubset.size(); j++) {
              u = nodeSubset[j];
              v = sibling[u];
              if (is_pickup[u] && !inS[v] && nodeSubset[j] != loaded_in_k[0]) {
                loaded_in_k.push_back(u);
              }
            }
            for (int kk = 0; kk < K; kk++) {
              ctValue1 = 0.0;
              IloRange cut = generate_lifted_lifo(paths[p][i], loaded_in_k[0], kk, ctValue1, expr, loaded_in_k);
              expr.clear();
              if (ctValue1 > nodeSubset.size() + 2 + EPSILON) {
                violatedLifos[paths[p][i]][loaded_in_k[0]] = true;
                lifos.add(cut);
                lifoCounter++;
                flag = true;
              }
            }
          }

          loaded_in_k.clear();
          nodeSubset.clear();
          memset(inS, false, sizeof (inS));
          possible_lifo = false;

          nodeSubset.push_back(paths[p][i]);
          loaded_in_k.push_back(paths[p][i]);

          int last_delivery = -1;
          for (int j = i + 1; j < last; ++j) {
            u = paths[p][j];
            v = sibling[u];
            nodeSubset.push_back(u);
            inS[u] = true;
            if (!is_pickup[u] && visited[u] && !visited[v]) {
              last_delivery = j;
              j = last;
              possible_lifo = true;
              nodeSubset.pop_back();
            }
          }

          if (possible_lifo && violatedLifos[sibling[paths[p][last_delivery]]][paths[p][i]]) {
            for (int j = 0; j < nodeSubset.size(); j++) {
              u = nodeSubset[j];
              v = sibling[u];
              if (is_pickup[u] && !inS[v] && nodeSubset[j] != loaded_in_k[0]) {
                loaded_in_k.push_back(u);
              }
            }
            for (int kk = 0; kk < K; kk++) {
              ctValue1 = 0.0;
              IloRange cut = generate_lifted_lifo(sibling[paths[p][last_delivery]], paths[p][i], kk, ctValue1, expr, loaded_in_k);
              expr.clear();
              if (ctValue1 > nodeSubset.size() + 2 + EPSILON) {
                violatedLifos[sibling[paths[p][last_delivery]]][paths[p][i]] = true;
                lifos.add(cut);
                lifoCounter++;
                flag = true;
              }
            }
          }

          loaded_in_k.clear();
          nodeSubset.clear();

        }


        if (original_lifo) {

          for (int j = last - 1; j >= i + 1; j--) {
            u = paths[p][j];
            v = sibling[u];
            nodeSubset.push_back(u);
            if (is_pickup[u]) {
              if (visited[u] && !visited[v]) {
                for (int kk = 0; kk < K; kk++) {
                  ctValue1 = 0.0;
                  IloRange cut = generateLIFOS(paths[p][i], u, kk, ctValue1, expr);
                  expr.clear();
                  if (ctValue1 > nodeSubset.size() + 2 + EPSILON) {
                    lifos.add(cut);
                    lifoCounter++;
                    j = 0;
                    flag = true;
                  }
                }
              }
            }
          }

          nodeSubset.clear();

          nodeSubset.push_back(paths[p][i]);

          for (int j = i + 1; j < last; j++) {
            u = paths[p][j];
            v = sibling[u];
            nodeSubset.push_back(u);
            if (!is_pickup[u]) {
              if (visited[u] && !visited[v]) {
                visited[u] = false;
                nodeSubset.pop_back();
                for (int kk = 0; kk < K; kk++) {
                  ctValue1 = 0.0;
                  IloRange cut = generateLIFOS(v, paths[p][i], kk, ctValue1, expr);
                  expr.clear();
                  if (ctValue1 > nodeSubset.size() + 2 + EPSILON) {
                    lifos.add(cut);
                    lifoCounter++;
                    j = last;
                    flag = true;
                  }
                }
              }
            }
          }

          nodeSubset.clear();

        }

        if (!flag) {//check if the set defines a violated concflict capacity inequality (37)
          nodeSubset.clear();
          for (int j = i + 1; j < last; j++)
            nodeSubset.push_back(paths[p][j]);
          /****************************************
           *******Conflict Capacity Constraints*****
           *****************************************/
          short int phi = 0;
          short int the = 0;
          short int zS = 0;
          vector<int> predecessors;
          vector<int> successors;
          memset(visited, false, sizeof (visited));
          double qs = 0;
          for (int l = 0; l < nodeSubset.size(); l++) {
            if (!is_pickup[nodeSubset[l]]) {
              predecessors.push_back(sibling[nodeSubset[l]]);
              qs -= 1; //needless
            } else {
              successors.push_back(sibling[nodeSubset[l]]);
              qs += 1; //needless
            }
            visited[nodeSubset[l]] = true;
          }
          for (int l = 0; l < predecessors.size(); l++)
            if (!visited[predecessors[l]])
              phi += itemLength[predecessors[l]];
          for (int l = 0; l < successors.size(); l++)
            if (!visited[successors[l]])
              the += -itemLength[successors[l]];
          zS = phi >= the ? phi : the;
          if (zS > (Q)*(K - 1) + EPSILON) {
            ctValue1 = 0.0;
            IloRange capCut = generateConflicitCapacity(paths[p][i], ctValue1, expr);
            caps.add(capCut);
            expr.clear();
            iCapCounter++;
          }
          nodeSubset.clear();
        }
      }
      u = paths[p][i];
    }

    u = 0;
    expr.clear();
    memset(capacities, 0, sizeof (capacities));
    memset(visited2, false, sizeof (visited2));
    nodeSubset.clear();

    if (sec_capacity) {

      for (int i = 1; i <= pIndex; i++) {
        v = paths[p][i];
        int k = used_stack[v];
        if (is_pickup[v]) {
          old_capacity = capacities[k];
          capacities[k] += itemLength[v];
          visited2[v] = true;
          if (old_capacity == 0) {//before loading v, stack k was empty...
            start[k] = i - 1;
          } else if (capacities[k] > Q) {
            end[k] = i;

            nodeSubset.push_back(paths[p][start[k]]);
            for (int ti = start[k]; ti < end[k]; ti++) {
              if (is_pickup[paths[p][ti + 1]] && used_stack[paths[p][ti + 1]] == k) {
                nodeSubset.push_back(paths[p][ti + 1]);
                loaded_in_k.push_back(paths[p][ti + 1]);
              } else {
                nodeSubset.push_back(paths[p][ti + 1]);
              }
            }

            for (int kk = 0; kk < K; kk++) {
              ctValue1 = 0.0;
              IloRange cut = generateSECCapacity(kk, ctValue1, expr, loaded_in_k);
              expr.clear();
              if (ctValue1 > loaded_in_k.size() + nodeSubset.size() - 2 + EPSILON) {
                caps.add(cut);
                pCapCounter++;
              }
            }

            loaded_in_k.clear();
            nodeSubset.clear();

          }
        }
        if (!is_pickup[v] && visited2[sibling[v]]) {
          capacities[k] -= -itemLength[v];
        }
      }

    }

    if (path_capacity) {

      int count = 0;
      ctValue1 = 0;
      for (int i = 1; i < pIndex; i++) {
        v = paths[p][i];
        int k = used_stack[v];
        if (is_pickup[v]) {
          old_capacity = capacities[k];
          capacities[k] += itemLength[v];
          visited2[v] = true;
          if (old_capacity == 0) //before loading v, stack k was empty...
          {
            start[k] = i - 1;
          } else if (capacities[k] > Q) {
            end[k] = i;
            for (int kk = 0; kk < K; kk++) {//generate a paths[p] inequality
              for (int ti = start[k]; ti < end[k]; ti++) {
                if (is_pickup[paths[p][ti + 1]] && used_stack[paths[p][ti + 1]] == k) {
                  expr += y[ mapVarsY[10 * paths[p][ti + 1] + kk] ];
                  expr += x[ mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ];
                  ctValue1 += _x[ mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ] + _y[ mapVarsY[10 * paths[p][ti + 1] + kk] ];
                  count += 2;
                } else {
                  count++;
                  expr += x[ mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ];
                  ctValue1 += _x[ mapVarsX[100 * paths[p][ti] + paths[p][ti + 1]] ];
                }
              }
              if (ctValue1 > count - 1 + EPSILON) {
                caps.add(expr <= count - 1);
                pCapCounter++;
              }
              expr.clear();
              count = 0;
              ctValue1 = 0;
            }
          }
        }
        if (!is_pickup[v] && visited2[sibling[v]]) {
          capacities[k] -= -itemLength[v];
        }
      }

    }

  }

#endif

  if (secs.getSize() == 0 && precs.getSize() == 0 && lifos.getSize() == 0 && caps.getSize() == 0) {
    flag = 0;
  } else {
    flag = true;
    counter1 = secs.getSize();
    counter2 = precs.getSize();
    counter3 = lifos.getSize();
    counter4 = caps.getSize();
  }

  return flag;
}

void TSPPDMS::addCuts() {
  try {
    if (secs.getSize() > 0) {
      Model.add(secs);
    }
    if (precs.getSize() > 0) {
      Model.add(precs);
    }
    if (lifos.getSize() > 0) {
      Model.add(lifos);
    }
    if (caps.getSize() > 0) {
      Model.add(caps);
    }
  } catch (IloException& ex) {
    cout << "addCuts: " << ex << endl;
  }
}

void TSPPDMS::modifyGraph(Dinic& graph, int pi, int dj, int stk) {
  int di = sibling[pi];
  int pj = sibling[dj];
  int kap, aux;
  bool flag;
  kap = aux = 0;
  flag = false;

  for (int u = 1; u < 2 * n; u++) { // agregar em (u,i)
    if (u != pi && u != pj) {
      if (_x[ mapVarsX[100 * u + pi] ] > 0.5) {
        aux = (int) (10000 * _x[ mapVarsX[100 * u + pi] ] + 0.5);
        for (int k = 0; k < K; k++) {
          if (k != stk && !(u == 0 && k != 0)) {
            kap = (int) (10000 * _y[ mapVarsY[10 * pi + k] ] + 0.5);
          }
        }
        if (kap > 100) {
          flag = graph.modifyArc(pi, u, kap + aux);
          if (!flag)
            graph.addArc(pi, u, kap);
        }
      }
    }
  }

  for (int u = 1; u < 2 * n; u++) { // agregar em (u,j)
    if (u != pj) {
      if (_x[ mapVarsX[100 * u + pj] ] > 0.5) {
        aux = (int) (10000 * _x[ mapVarsX[100 * u + pj] ] + 0.5);
        for (int k = 0; k < K; k++) {
          if (k != stk && !(u == 0 && k != 0)) {
            kap = (int) (10000 * _y[ mapVarsY[10 * pj + k] ] + 0.5);
          }
        }
        if (kap > 100) {
          flag = graph.modifyArc(u, pj, kap + aux);
          if (!flag)
            graph.addArc(u, pj, kap);
        }
      }
    }
  }

  for (int u = 1; u < 2 * n; u++) { // agregar em (u,j)
    if (u != dj) {
      if (_x[ mapVarsX[100 * u + dj] ] > 0.5) {
        flag = graph.modifyArc(pi, u, 10000);
        if (!flag) {
          graph.addArc(pi, u, 10000);
        }
      }
    }
  }

  graph.addArc(dj, pj, 50000);
  flag = graph.modifyArc(di, pj, 50000);
  if (!flag) {
    graph.addArc(di, pj, 50000);
  }

}

int TSPPDMS::setModelProblem() {
  int idx;
  Arc* e;
  char name[25];
  int posy, posx;
  int pi, pj, di, dj;

  y = IloNumVarArray(Env, K * n, 0.0, 1.0, ILOINT);
  _y = IloNumArray(Env, K * n);
  x = IloNumVarArray(Env, 4 * n * n - n, 0.0, 1.0, ILOINT);
  _x = IloNumArray(Env, 4 * n * n - n);

  convy = IloConversion(Env, y, ILOFLOAT);
  convx = IloConversion(Env, x, ILOFLOAT);

  Constraints = IloRangeArray(Env);

  /* Objective Funtion */
  IloExpr xfo(Env);
  posy = posx = 0;

  for (int i = 0; i < V; i++) {
    for (int ii = 0; ii < adj[i].size(); ii++) {
      e = adj[i].at(ii);
      idx = e->u != i ? e->u : e->v;
      if (!(!is_pickup[i] && sibling[i] == idx) && idx != 0 && i != 2 * n + 1
              && !(is_pickup[i] && idx == 2 * n + 1) && !(i == 0 && !is_pickup[idx])) {
        sprintf(name, "x(%d,%d)", i, idx);
        x[posx].setName(name);
        xfo += (e->c) * x[posx];
        mapVarsX[100 * i + idx] = posx++;
      }
    }
  }

  fo = IloAdd(Model, IloMinimize(Env, xfo));
  xfo.end();

  /* Y variables*/
  for (int i = 0; i < pickups.size(); ++i) {
    for (int k = 0; k < K; ++k) {
      sprintf(name, "y(%d,%d)", pickups[i], k);
      y[posy].setName(name);
      mapVarsY[10 * pickups[i] + k] = posy++;
    }
  }

  IloExpr expr(Env);

  /* Degree constraints */

  for (int i = 1; i <= 2 * n; i++) {
    for (int j = 1; j <= 2 * n; j++) {
      if (i != j && !(!is_pickup[i] && sibling[i] == j)) {
        expr += x[ mapVarsX[100 * i + j] ];
      }
    }
    if (!is_pickup[i])
      expr += x[ mapVarsX[100 * i + (2 * n + 1)] ];
    sprintf(name, "degOut_%d", i);
    IloRange crt;
    crt = 1 <= expr <= 1;
    crt.setName(name);
    Constraints.add(crt);
    expr.clear();
  }
  for (int j = 1; j <= 2 * n; j++) {
    for (int i = 1; i <= 2 * n; i++) {
      if (i != j && !(!is_pickup[i] && sibling[i] == j)) {
        expr += x[ mapVarsX[100 * i + j] ];
      }
    }
    if (is_pickup[j])
      expr += x[ mapVarsX[100 * 0 + j] ];
    sprintf(name, "degIn_%d", j);
    IloRange crt;
    crt = 1 <= expr <= 1;
    crt.setName(name);
    Constraints.add(crt);
    expr.clear();
  }


  //degree constraint for depot 0
  for (int i = 0; i < pickups.size(); i++)
    expr += x[ mapVarsX[pickups[i]] ];
  sprintf(name, "degOut_0");
  IloRange crt1;
  crt1 = 1 <= expr <= 1;
  crt1.setName(name);
  Constraints.add(crt1);
  expr.clear();

  //degree constraint for depot 2n+1
  for (int i = 0; i < deliveries.size(); i++)
    expr += x[ mapVarsX[100 * deliveries[i] + 2 * n + 1] ];
  sprintf(name, "degIn_%d", 2 * n + 1);
  IloRange crt2;
  crt2 = 1 <= expr <= 1;
  crt2.setName(name);
  Constraints.add(crt2);
  expr.clear();

  /*The successor of depot 0 is a pickup on stack 0*/
  for (int i = 0; i < pickups.size(); ++i) {
    expr = x[mapVarsX[100 * 0 + pickups[i]]] - y[mapVarsY[10 * pickups[i] + 0]];
    IloRange crt;
    crt = expr <= 0;
    sprintf(name, "sym_%d", pickups[i]);
    crt.setName(name);
    Constraints.add(crt);
    expr.clear();
  }

  /* Stack constraints*/

  for (int i = 0; i < pickups.size(); i++) {
    for (int k = 0; k < K; ++k) {
      expr += y[ mapVarsY[10 * pickups[i] + k] ];
    }
    sprintf(name, "stk_%d", pickups[i]);
    IloRange crt;
    crt = (1 <= expr <= 1);
    crt.setName(name);
    Constraints.add(crt);
    expr.clear();
  }

  expr.clear();

  /* SECS with |S|=2 */
  for (int i = 1; i <= 2 * n; i++) {
    for (int j = i + 1; j <= 2 * n; j++) {
      if (!(!is_pickup[i] && sibling[i] == j) && !(!is_pickup[j] && sibling[j] == i)
              && !(is_pickup[i] && !is_pickup[j]) && !(!is_pickup[i] && is_pickup[j])) {
        expr += x[ mapVarsX[100 * i + j] ] + x[ mapVarsX[100 * j + i] ];
        sprintf(name, "sec2_%d%d", i, j);
        IloRange sec2;
        sec2 = expr <= 1;
        sec2.setName(name);
        Constraints.add(sec2);
        expr.clear();
      }
    }
  }

  for (int i = 0; i < pickups.size(); i++) {
    for (int j = i + 1; j < pickups.size(); j++) {
      int pi = pickups[i];
      int pj = pickups[j];
      int di = sibling[pi];
      int dj = sibling[pj];
      if (di != dj) {
        expr += x[ mapVarsX[100 * pi + dj] ] + x[ mapVarsX[100 * dj + pi] ] + x[ mapVarsX[100 * pj + di] ] + x[ mapVarsX[100 * di + pj] ];
        Constraints.add(expr <= 1);
        expr.clear();
      }
    }
  }

#ifdef SUCC_INEQ
  for (int i = 0; i < pickups.size(); i++) {
    int pi = pickups[i];
    int di = sibling[pi];
    for (int j = 0; j < pickups.size(); j++) {
      int pj = pickups[j];
      int dj = sibling[pj];
      if (pi == pj)
        continue;
      for (int k = 0; k < 1; k++) {
        expr = y[ mapVarsY[10 * pi + k] ] + y[ mapVarsY[10 * pj + k] ] + x[ mapVarsX[ 100 * pi + dj ] ];
        sprintf(name, "succ_ineq_%d_%d", pi, pj);
        IloRange succ;
        succ = expr <= 2;
        succ.setName(name);
        Constraints.add(succ);
        expr.clear();
      }
    }
  }
#endif

#ifdef INITIAL_PRECEDENCE_S1

  // Only one pickup in S

  for (int i = 0; i < pickups.size(); ++i) {
    expr = x[mapVarsX[100 * 0 + pickups[i]]];
    for (int j = 0; j < deliveries.size(); ++j) {
      if (deliveries[j] != sibling[pickups[i]]) {
        expr += x[mapVarsX[100 * pickups[i] + deliveries[j]]];
      }
    }
    IloRange ilorng;
    sprintf(name, "s_initial_prec_%d", pickups[i]);
    ilorng = expr <= 1;
    Constraints.add(ilorng);
    expr.clear();
  }

  for (int i = 0; i < deliveries.size(); ++i) {
    expr = x[mapVarsX[100 * deliveries[i] + (2 * n + 1)]];
    for (int j = 0; j < pickups.size(); ++j) {
      if (deliveries[i] != sibling[pickups[j]]) {
        expr += x[mapVarsX[100 * pickups[j] + deliveries[i]]];
      }
    }
    IloRange ilorng;
    sprintf(name, "e_initial_prec_%d", deliveries[i]);
    ilorng = expr <= 1;
    ilorng.setName(name);
    Constraints.add(ilorng);
    expr.clear();
  }

#endif  

#ifdef INITIAL_PRECEDENCE_S2

  // Two pickups in S

  // Route Start Inequalities
  for (int i = 0; i < pickups.size(); ++i) {
    pi = pickups[i];
    di = sibling[pi];
    for (int j = j + 1; j < pickups.size(); j++) {
      pj = pickups[j];
      dj = sibling[pj];
      if (pi != pj) {
        expr = x[mapVarsX[100 * 0 + pi]] + x[mapVarsX[100 * 0 + pj]];
        expr += x[mapVarsX[100 * pi + pj]] + x[mapVarsX[100 * pj + pi]];
        for (int l = 0; l < deliveries.size(); ++l) {
          if (deliveries[l] != di && deliveries[l] != dj) {
            expr += x[mapVarsX[100 * pi + deliveries[l]]];
            expr += x[mapVarsX[100 * pj + deliveries[l]]];
          }
        }
        IloRange ilorng;
        sprintf(name, "s_initial_prec_%d_%d", pi, pj);
        ilorng = expr <= 2;
        ilorng.setName(name);
        Constraints.add(ilorng);
        expr.clear();
      }
    }
  }

  // Route End Inequalities 
  for (int i = 0; i < pickups.size(); ++i) {
    pi = pickups[i];
    di = sibling[pi];
    for (int j = 0; j < pickups.size(); j++) {
      pj = pickups[j];
      dj = sibling[pj];
      if (pi != pj) {
        expr = x[mapVarsX[100 * di + (2 * n + 1)]] + x[mapVarsX[100 * dj + (2 * n + 1)]];
        expr += x[mapVarsX[100 * di + dj]] + x[mapVarsX[100 * dj + di]];
        for (int l = 0; l < pickups.size(); ++l) {
          if (pickups[l] != pi && pickups[l] != pj) {
            expr += x[mapVarsX[100 * pickups[l] + di]];
            expr += x[mapVarsX[100 * pickups[l] + dj]];
          }
        }
        IloRange ilorng;
        sprintf(name, "e_initial_prec_%d_%d", di, dj);
        ilorng = expr <= 2;
        ilorng.setName(name);
        Constraints.add(ilorng);
        expr.clear();
      }
    }
  }

#endif  

#ifdef INITIAL_PRECEDENCE_S3

  // Three pickups in S

  int pl, dl;

  // Route Start Inequalities
  for (int i = 0; i < pickups.size(); ++i) {
    pi = pickups[i];
    di = sibling[pi];
    for (int j = i + 1; j < pickups.size(); j++) {
      pj = pickups[j];
      dj = sibling[pj];
      if (pi != pj) {
        for (int l = j + 1; l < pickups.size(); l++) {
          pl = pickups[l];
          dl = sibling[pl];
          if (pi != pl && pl != pj) {
            expr = x[mapVarsX[100 * 0 + pi]] + x[mapVarsX[100 * 0 + pj]] + x[mapVarsX[100 * 0 + pl]];
            expr += x[mapVarsX[100 * pi + pj]] + x[mapVarsX[100 * pj + pl]];
            expr += x[mapVarsX[100 * pi + pl]] + x[mapVarsX[100 * pl + pj]];
            expr += x[mapVarsX[100 * pj + pi]];
            expr += x[mapVarsX[100 * pl + pi]];
            for (int a = 0; a < deliveries.size(); ++a) {
              if (deliveries[a] != di && deliveries[a] != dj && deliveries[a] != dl) {
                expr += x[mapVarsX[100 * pi + deliveries[a]]];
                expr += x[mapVarsX[100 * pj + deliveries[a]]];
                expr += x[mapVarsX[100 * pl + deliveries[a]]];
              }
            }
            IloRange ilorng;
            sprintf(name, "initial_prec_%d_%d", pi, pj);
            ilorng = expr <= 3;
            Constraints.add(ilorng);
            expr.clear();
          }
        }
      }
    }
  }

  // Route End Inequalities 
  for (int i = 0; i < pickups.size(); ++i) {
    pi = pickups[i];
    di = sibling[pi];
    for (int j = i + 1; j < pickups.size(); j++) {
      pj = pickups[j];
      dj = sibling[pj];
      if (pi != pj) {
        for (int l = j + 1; l < pickups.size(); l++) {
          pl = pickups[l];
          dl = sibling[pl];
          if (pl != pi && pl != pj) {
            expr = x[mapVarsX[100 * di + (2 * n + 1)]] + x[mapVarsX[100 * dj + (2 * n + 1)]] + x[mapVarsX[100 * dl + (2 * n + 1)]];
            expr += x[mapVarsX[100 * di + dj]] + x[mapVarsX[100 * di + dl]];
            expr += x[mapVarsX[100 * dj + di]] + x[mapVarsX[100 * dj + dl]];
            expr += x[mapVarsX[100 * dl + di]] + x[mapVarsX[100 * dl + dj]];
            for (int a = 0; a < pickups.size(); ++a) {
              if (pickups[a] != pi && pickups[a] != pj && pickups[a] != pl) {
                expr += x[mapVarsX[100 * pickups[a] + di]];
                expr += x[mapVarsX[100 * pickups[a] + dj]];
                expr += x[mapVarsX[100 * pickups[a] + dl]];
              }
            }
            IloRange ilorng;
            sprintf(name, "initial_prec_%d_%d", di, dj);
            ilorng = expr <= 3;
            Constraints.add(ilorng);
            expr.clear();
          }
        }
      }
    }

  }

#endif   


  Model.add(Constraints);
  //Cplex.exportModel("initial_model.lp");

  secs = IloRangeArray(Env);
  precs = IloRangeArray(Env);
  lifos = IloRangeArray(Env);
  caps = IloRangeArray(Env);

  IloNumArray priorities(Env, 4 * n * n - n);
  for (int i = 0; i < 4 * n * n - n; i++)
    priorities[i] = 2.0;
  Cplex.setPriorities(x, priorities);

  //  cout << "Success setting model." << endl;

  return 0;
}

int TSPPDMS::solveProblem() {
  try {
    Cplex.solve();
    Cplex.getValues(_x, x);
    Cplex.getValues(_y, y);
    if (rlModel)
      root = Cplex.getObjValue();
  } catch (IloException& ex) {
    cout << "solveProblem:" << ex << endl;
  }
  return 0;
}

void TSPPDMS::relaxIntVars() {
  Model.add(convy);
  Model.add(convx);
  rlModel = true;
  return;
}

void TSPPDMS::enforceIntVars() {
  Model.remove(convy); //integral x.
  Model.remove(convx); //integral y.
  rlModel = false;
  //Cplex.exportModel("testRL.lp");
  return;
}

IloNum TSPPDMS::getSolution() {
  return Cplex.getObjValue();
}

void TSPPDMS::addCut(IloRange cut) {
  cout << cut << endl << endl;
}

int TSPPDMS::getNumRequisition() {
  return n;
}

/*
 x(S) <= |S| - 1
 */
IloRange TSPPDMS::generateSEC(IloNum& val, IloExpr & cut) {
  int i, j;
  int m;
  int idx;

  m = nodeSubset.size();

  for (int ii = 0; ii < m; ii++)
    for (int jj = 0; jj < m; jj++) {
      i = nodeSubset[ii];
      j = nodeSubset[jj];
      if (i != j && !(!is_pickup[i] && sibling[i] == j) && j != 0 && i != 2 * n + 1
              && !(is_pickup[i] && j == 2 * n + 1) && !(i == 0 && !is_pickup[j])) {
        idx = mapVarsX[100 * i + j];
        cut += x[idx];
        val += _x[idx];
      }
    }

  m--;


#ifdef LIFTED_PRECEDENCE_SEC
  bool incS[2 * n + 2];
  memset(incS, false, sizeof (incS));
  for (int ii = 0; ii < nodeSubset.size(); ii++)
    incS[nodeSubset[ii]] = true;

  int pi, di, pj, dj;
  pi = di = pj = dj = -1;

  vector<int> d_of_lone_p;

  for (int i = 0; i < nodeSubset.size(); i++) {
    if (is_pickup[nodeSubset[i]] && !incS[sibling[nodeSubset[i]]]) {
      d_of_lone_p.push_back(sibling[nodeSubset[i]]);
    }
  }

  for (int i = 0; i < nodeSubset.size(); i++) {
    if (!is_pickup[nodeSubset[i]] && !incS[sibling[nodeSubset[i]]]) {
      pj = sibling[nodeSubset[i]];
      i = nodeSubset.size();
    }
  }

  if (d_of_lone_p.size() > 0 && pj != -1) {
    for (int i = 0; i < d_of_lone_p.size(); i++) {
      cut += x[mapVarsX[100 * pj + d_of_lone_p[i]]]; //x[mapVarsX[100 * d_of_lone_p[i] + pj]] 
    }
  }

#endif


  /*printf("##########################\n");
  printf("S = { ");
  for (int i = 0; i < nodeSubset.size(); i++) {
          printf("%d ", nodeSubset[i]);
  }
  printf("}\n");
  cout << cut << " <= " << m << endl;
  printf("##########################\n");*/

  return (cut <= m);
}

/*
 x(S) <= |S| - 2
 */
IloRange TSPPDMS::generatePREC(IloNum& val, IloExpr & cut) {
  int i, j;
  int m;
  int idx;

  m = nodeSubset.size();

  for (int ii = 0; ii < m; ii++)
    for (int jj = 0; jj < m; jj++) {
      i = nodeSubset[ii];
      j = nodeSubset[jj];
      if (i != j && !(!is_pickup[i] && sibling[i] == j) && j != 0 && i != 2 * n + 1
              && !(is_pickup[i] && j == 2 * n + 1) && !(i == 0 && !is_pickup[j])) {
        idx = mapVarsX[100 * i + j];
        cut += x[idx];
        val += _x[idx];
      }
    }

  m -= 2;
  return (cut <= m);
}

/*
 x(n+i,S) + x(S) + x(S,i) <= |S| : \forall S \subseteq V \setminus \{0,2n+1,i,n+i\}
 */
IloRange TSPPDMS::generatePRECI(IloNum& val, IloExpr& cut, int pi) {
  int i, j;
  int m;
  int di;

  di = sibling[pi];
  m = nodeSubset.size();

  //i,n+i,0,2n+1 \notin S
  //x(n+i,S)
  for (int t = 0; t < m; t++) {
    j = nodeSubset[t];
    cut += x[ mapVarsX[100 * di + j] ];
    val += _x[ mapVarsX[100 * di + j] ];
  }

  //x(S)
  for (int t = 0; t < m; t++) //in S
  {
    for (int u = 0; u < m; u++) //in S
    {
      if (t != u) {
        i = nodeSubset[t];
        j = nodeSubset[u];
        if (is_pickup[i] || sibling[i] != j)//exclude x_{n+i,i,k} from model)
        {
          cut += x[ mapVarsX[100 * i + j] ];
          val += _x[ mapVarsX[100 * i + j] ];
        }
      }
    }
  }

  //x(S,i)
  for (int t = 0; t < m; t++) {
    j = nodeSubset[t];
    cut += x[ mapVarsX[100 * j + pi] ];
    val += _x[ mapVarsX[100 * j + pi] ];
  }

  return (cut <= m);
}

/*
x(S) + x(S,\overline{S} \cap \pi(S)) + x(S \cap \pi(S),\overline{S} \setminus \pi(S)) <= |S| - 1
 */
IloRange TSPPDMS::generatePredecessorCut(IloNum& val, IloExpr & cut) {
  vector<int> pre;
  bool inS[2 * n + 2]; //S incidence vector
  bool inP[2 * n + 2]; //Pre incidence vector
  int s, t;
  int m;

  for (int i = 0; i <= 2 * n + 1; i++) {
    inS[i] = false;
    inP[i] = false;
  }

  m = nodeSubset.size();
  for (int i = 0; i < m; i++) {
    s = nodeSubset[i];
    inS[s] = true;
    if (!is_pickup[s]) //s is a delivery
    {
      pre.push_back(sibling[s]); //predecessor set: i\in P; n+i \in S
      inP[sibling[s]] = true;
    }

    for (int j = 0; j < m; j++) {
      if (i != j) {
        t = nodeSubset[j];
        if (!(!is_pickup[s] && sibling[s] == t) && t != 0 && s != 2 * n + 1
                && !(is_pickup[s] && t == 2 * n + 1) && !(s == 0 && !is_pickup[t])) {
          cut += x[ mapVarsX[100 * s + t] ];
          val += _x[ mapVarsX[100 * s + t] ];
        }
      }
    }
  }

  for (int i = 0; i < m; i++) //x_ij; i \in S j \in S'^Pre(S)
  {
    s = nodeSubset[i];
    for (int j = 0; j < pre.size(); j++) {
      t = pre[j];
      if (!inS[t] && !(!is_pickup[s] && sibling[s] == t) && t != 0 && s != 2 * n + 1
              && !(is_pickup[s] && t == 2 * n + 1) && !(s == 0 && !is_pickup[t])) {
        cut += x[ mapVarsX[100 * s + t] ];
        val += _x[ mapVarsX[100 * s + t] ];
      }
    }
  }

  for (int i = 0; i < m; i++) //x_ij; i in S^pre(S) j \in S'-pre(S)
  {
    s = nodeSubset[i];
    if (inP[s]) {
      for (int j = 1; j <= 2 * n + 1; j++)
        if ((!inS[j] && !inP[j]) && !(!is_pickup[s] && sibling[s] == j) && j != 0
                && s != 2 * n + 1 && !(is_pickup[s] && j == 2 * n + 1) && !(s == 0 && !is_pickup[j])) {
          cut += x[ mapVarsX[100 * s + j] ];
          val += _x[ mapVarsX[100 * s + j] ];
        }
    }
  }

  m = m - 1;
  return (cut <= m);
}

/*
 x(S) + x(\overline{S} \cap \sigma(S),S) + x(\overline{S} \setminus \sigma(S), S \cap \sigma(S)) <= |S| - 1
 */
IloRange TSPPDMS::generateSuccessorCut(IloNum& val, IloExpr & cut) {
  vector<int> succ;
  bool incS[2 * n + 2]; //S incidence vector
  bool incC[2 * n + 2]; //Suc incidence vector
  int s, t;
  int m;

  for (int i = 0; i <= 2 * n + 1; i++) {
    incS[i] = false;
    incC[i] = false;
  }

  m = nodeSubset.size();
  for (int i = 0; i < m; i++) //x(S)
  {
    s = nodeSubset[i];
    incS[s] = true;
    if (is_pickup[s]) //s is a pickup
    {
      succ.push_back(sibling[s]); //successor set: n+i\in D; i \in S
      incC[sibling[s]] = true;
    }

    for (int j = 0; j < m; j++) {
      if (i != j) {
        t = nodeSubset[j];
        if ((is_pickup[s] || sibling[s] != t) && t != 0 && s != 2 * n + 1
                && (!is_pickup[s] || t != 2 * n + 1) && (s != 0 || is_pickup[t])) {
          cut += x[ mapVarsX[100 * s + t] ];
          val += _x[ mapVarsX[100 * s + t] ];
        }
      }
    }
  }

  for (int j = 0; j < succ.size(); j++) {
    s = succ[j];
    if (!incS[s])
      for (int i = 0; i < m; i++) {
        t = nodeSubset[i];
        if ((is_pickup[s] || sibling[s] != t) && t != 0 && s != 2 * n + 1
                && (!is_pickup[s] || t != 2 * n + 1) && (s != 0 || is_pickup[t])) {
          cut += x[ mapVarsX[100 * s + t] ];
          val += _x[ mapVarsX[100 * s + t] ];
        }
      }
  }

  for (int j = 0; j < 2 * n + 1; j++) {
    if (!incS[j] && !incC[j])
      for (int i = 0; i < succ.size(); i++) {
        t = succ[i];
        if (incS[t] && (is_pickup[j] || sibling[j] != t) && t != 0 && j != 2 * n + 1
                && (!is_pickup[j] || t != 2 * n + 1) && (j != 0 || is_pickup[t])) {
          cut += x[ mapVarsX[100 * j + t] ];
          val += _x[ mapVarsX[100 * j + t] ];
        }
      }
  }

  m = m - 1;
  return (cut <= m);
}

IloRange TSPPDMS::generate_lifted_lifo(int pi, int pj, int k, IloNum &val, IloExpr &cut, vector<int> &loaded_in_k) {
  // First element in loaded_in_k is the first element visited in S
  int di = sibling[pi];
  int dj = sibling[loaded_in_k[0]];
  int m = nodeSubset.size();
  int i, j;
  int idx;
  bool incS[2 * n + 2];

  memset(incS, false, sizeof (incS));
  for (int ii = 0; ii < m; ii++)
    incS[nodeSubset[ii]] = true;

  idx = mapVarsY[10 * pi + k]; // x(\overline{S},loaded_in_k[0])
  val += _y[idx];
  cut += y[idx]; // y_{i}_^{k}
  for (int i = 0; i < loaded_in_k.size(); ++i) {
    idx = mapVarsY[10 * loaded_in_k[i] + k]; // x(\overline{S},loaded_in_k[0])
    val += _y[idx] / loaded_in_k.size();
    cut += y[idx] / loaded_in_k.size(); // y_{j}^{k} : j \in loaded_in_k
  }

  for (i = 1; i < 2 * n + 1; i++) {
    if (!incS[i] && i != dj && i != loaded_in_k[0] && i != di) {
      idx = mapVarsX[100 * i + loaded_in_k[0]]; // x(\overline{S},loaded_in_k[0])
      cut += x[idx];
      val += _x[idx];
    }
  }

  for (int ii = 0; ii < m; ii++) { // x(S)
    for (int jj = 0; jj < m; jj++) {
      i = nodeSubset[ii];
      j = nodeSubset[jj];
      if (i != j && (is_pickup[i] || sibling[i] != j) && j != 0 && i != 2 * n + 1 && (!is_pickup[i] || j != 2 * n + 1) && (i != 0 || is_pickup[j])) {
        idx = mapVarsX[100 * i + j];
        cut += x[idx];
        val += _x[idx];
      }
    }
  }

  for (int t = 0; t < m; t++) {
    idx = mapVarsX[100 * nodeSubset[t] + di];
    cut += x[idx];
    val += _x[idx];
  }
  return (cut <= m + 2);
}

IloRange TSPPDMS::generate_lifted_lifo_loading_inside(int pi, int pj, int k, IloNum &val, IloExpr & cut) {
  // First element in loaded_in_k is the first element visited in S


  int di = sibling[pi];
  int dj = sibling[pj];
  int m = nodeSubset.size();
  int i, j;
  int idx;
  bool incS[2 * n + 2];

  vector<int> loaded_in_k;

  memset(incS, false, sizeof (incS));
  for (int ii = 0; ii < m; ii++)
    incS[nodeSubset[ii]] = true;

  int u, v;
  for (int j = 0; j < m; ++j) {
    u = nodeSubset[j];
    v = sibling[u];
    if (is_pickup[u] && incS[u] && !incS[v]) {
      loaded_in_k.push_back(u);
    }
  }

  idx = mapVarsY[10 * pi + k]; // x(\overline{S},loaded_in_k[0])
  val += _y[idx];
  cut += y[idx]; // y_{i}_^{k}
  for (int i = 0; i < loaded_in_k.size(); ++i) {
    idx = mapVarsY[10 * loaded_in_k[i] + k]; // x(\overline{S},loaded_in_k[0])
    val += _y[idx] / loaded_in_k.size();
    cut += y[idx] / loaded_in_k.size(); // y_{j}^{k} : j \in loaded_in_k
  }

  for (i = 1; i < 2 * n + 1; i++) {
    if (!incS[i] && i != dj && i != pj && i != di) {
      idx = mapVarsX[100 * i + pj]; // x(\overline{S},loaded_in_k[0])
      cut += x[idx];
      val += _x[idx];
    }
  }

  for (int ii = 0; ii < m; ii++) { // x(S)
    for (int jj = 0; jj < m; jj++) {
      i = nodeSubset[ii];
      j = nodeSubset[jj];
      if (i != j && (is_pickup[i] || sibling[i] != j) && j != 0 && i != 2 * n + 1 && (!is_pickup[i] || j != 2 * n + 1) && (i != 0 || is_pickup[j])) {
        idx = mapVarsX[100 * i + j];
        cut += x[idx];
        val += _x[idx];
      }
    }
  }

  for (int t = 0; t < m; t++) {
    idx = mapVarsX[100 * nodeSubset[t] + di];
    cut += x[idx];
    val += _x[idx];
  }
  return (cut <= m + 2);
}

/*
x(\overline{S},j) + x(S) + x(\overline{S},i) + x(S,n+j) <= |S| + 1 : j \in S, i, n+i, n+j \notin S
 */

IloRange TSPPDMS::generateLIFOS(int pi, int pj, int k, IloNum &val, IloExpr & cut) {
  int di = sibling[pi];
  int dj = sibling[pj];
  int m = nodeSubset.size();
  int i, j;
  int idx;
  bool incS[2 * n + 2];

  memset(incS, false, sizeof (incS));
  for (int ii = 0; ii < m; ii++)
    incS[nodeSubset[ii]] = true;

  for (i = 1; i < 2 * n + 1; i++) {
    if (!incS[i] && i != dj && i != pj && i != di) {
      idx = mapVarsX[100 * i + pj];
      cut += x[idx];
      val += _x[idx];
    }
  }

  idx = mapVarsY[10 * pj + k];
  cut += y[idx];
  val += _y[idx];

  for (int ii = 0; ii < m; ii++)
    for (int jj = 0; jj < m; jj++) {
      i = nodeSubset[ii];
      j = nodeSubset[jj];
      if (i != j && (is_pickup[i] || sibling[i] != j) && j != 0 && i != 2 * n + 1 && (!is_pickup[i] || j != 2 * n + 1) && (i != 0 || is_pickup[j])) {
        idx = mapVarsX[100 * i + j];
        cut += x[idx];
        val += _x[idx];
      }
    }

  idx = mapVarsY[10 * pi + k];
  cut += y[idx];
  val += _y[idx];

  for (int t = 0; t < m; t++) {
    pj = nodeSubset[t];
    idx = mapVarsX[100 * pj + di];
    cut += x[idx]; //di\notin S, pi \notin S x{di,i} do not appear here!
    val += _x[idx];
  }

  return (cut <= m + 2);
}

/*
 x(i,S) + x(S) + x(S,n+i) <= |S| : \forall S \subset P \cup D, i, n+i \notin S, z(S) > Q(K-1)
 */
IloRange TSPPDMS::generateConflicitCapacity(int pi, IloNum &val, IloExpr & cut) {
  int i, j;
  int m;
  int di;

  di = sibling[pi];
  m = nodeSubset.size();

  //{0,2n+1, pi, di}\notin S

  //y(i,S)
  for (int t = 0; t < m; t++) {
    j = nodeSubset[t];
    if (j != pi) {
      cut += x[ mapVarsX[100 * pi + j] ];
      val += _x[ mapVarsX[100 * pi + j] ];
    }
  }

  //y(S)
  for (int t = 0; t < m; t++) //in S
  {
    i = nodeSubset[t];
    for (int u = 0; u < m; u++) //in S
    {
      if (t != u) {
        j = nodeSubset[u];
        if (is_pickup[i] || sibling[i] != j)//exclude x_{n+i,i,k} from model)
        {
          cut += x[ mapVarsX[100 * i + j] ];
          val += _x[ mapVarsX[100 * i + j] ];
        }
      }
    }
  }

  //y(n+i,S)
  for (int t = 0; t < m; t++) {
    j = nodeSubset[t];
    if (j != di) {
      cut += x[ mapVarsX[100 * j + di] ];
      val += _x[ mapVarsX[100 * j + di] ];
    }
  }

  return (cut <= m);
}

/*
 x(S,\overline{S}) + x(\overline{S},S) >= \lceil |q(S)|/KQ \rceil
 */
IloRange TSPPDMS::generateRoundedCapacity(IloNum &val, IloExpr & cut) {
  int m;
  int lb;
  bool incS[MAXNODES];
  int pi, pj;

  //0,2n+1 \notin S
  memset(incS, false, sizeof (incS));

  m = nodeSubset.size();
  for (int i = 0; i < m; i++) {
    incS[nodeSubset[i]] = true;
  }

  double piS=0;
  double sigmaS=0;
  for (int i = 0; i < m; ++i) {
    if (!is_pickup[nodeSubset[i]] && !incS[sibling[nodeSubset[i]]]){
      piS += itemLength[sibling[nodeSubset[i]]];
    } else if (is_pickup[nodeSubset[i]] && !incS[sibling[nodeSubset[i]]]) {
      sigmaS += itemLength[nodeSubset[i]];
    }
  }

  double qs=max(piS,sigmaS);
  
  for (int i = 0; i < m; i++) {
    for (int j = 1; j <= 2 * n + 1; j++) {
      pi = nodeSubset[i];
      if (!incS[j] && pi != j && !(!is_pickup[pi] && sibling[pi] == j) && !(is_pickup[pi] && j == (2 * n + 1)))//delta_+
      {
        cut += x[ mapVarsX[100 * pi + j] ];
        val += _x[ mapVarsX[100 * pi + j] ];
      }
    }
  }

  /*for (int j = 0; j < m; j++) {
    for (int i = 0; i < 2 * n + 1; i++) {
      pj = nodeSubset[j];
      if (!incS[i] && pj != i && !(!is_pickup[i] && sibling[i] == pj) && !(i == 0 && !is_pickup[pj]))//delta_-
      {
        cut += x[ mapVarsX[100 * i + pj] ];
        val += _x[ mapVarsX[100 * i + pj] ];
      }
    }
  }*/

  lb = ceil(qs / (K * Q));

  //cout<<"lower bound:"<<qs<<"/"<<K*Q<<"="<<lb<<endl;

  return (cut >= lb);
}

IloRange TSPPDMS::generateSECCapacity(int k, IloNum &val, IloExpr &cut, vector<int> &loaded_in_k) {
  int m = nodeSubset.size();
  int t = loaded_in_k.size();
  int i, j;
  int idx;

  for (int ii = 0; ii < m; ii++) {
    for (int jj = 0; jj < m; jj++) {
      i = nodeSubset[ii];
      j = nodeSubset[jj];
      if (i != j && (is_pickup[i] || sibling[i] != j) && j != 0 && i != 2 * n + 1 && (!is_pickup[i] || j != 2 * n + 1) && (i != 0 || is_pickup[j])) {
        idx = mapVarsX[100 * i + j];
        cut += x[idx];
        val += _x[idx];
      }
    }
  }

  for (int l = 0; l < loaded_in_k.size(); l++) {
    idx = mapVarsY[10 * loaded_in_k[l] + k];
    cut += y[idx];
    val += _y[idx];
  }

  return (cut <= ((m - 1) + t) - 1); // m-1:path ; t:all items loaded in k
}


void TSPPDMS::setCplexSettings(int vlDisp, int vlEmph, int alg, int numThreads, double vlGap, int mem) {
  Cplex.use(lazyCallback(Env, *this)); //, _y, _x));       //lazy Constraints
  Cplex.use(cutCallback(Env, *this)); //, _y, _x));      //userDefined Constraints
  Cplex.use(incumbentCallback(Env, *this)); //, _y, _x));      //
  //Cplex.use(MipCallback(Env));

  Cplex.setParam(IloCplex::CutUp, initial_ub); //Sets the upper cutoff tolerance.

  if (vlDisp != 0) {
    Cplex.setParam(IloCplex::MIPDisplay, 2);
    Cplex.setParam(IloCplex::MIPInterval, vlDisp);
  } else
    Cplex.setOut(Env.getNullStream());
  Cplex.setWarning(Env.getNullStream());

  Cplex.setParam(IloCplex::RootAlg, alg);
  Cplex.setParam(IloCplex::NodeAlg, alg);
  Cplex.setParam(IloCplex::Threads, numThreads);
  Cplex.setParam(IloCplex::TiLim, 10800);
  if (vlEmph != 0)
    Cplex.setParam(IloCplex::MIPEmphasis, vlEmph);
  if (vlGap != 0.0) {
    Cplex.setParam(IloCplex::EpAGap, 1 - vlGap);
    Cplex.setParam(IloCplex::EpInt, vlGap);
    Cplex.setParam(IloCplex::ObjDif, 1 - vlGap);
    //Cplex.setParam(IloCplex::EpOpt,vlGap);
    //Cplex.setParam(IloCplex::EpGap,vlGap);
  }

  if (mem != 0)
    Cplex.setParam(IloCplex::WorkMem, mem);

  Cplex.setParam(IloCplex::MemoryEmphasis, 1);

  //Cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, 3); // Strong branching
  //Cplex.setParam(IloCplex::VarSel, 3);
  //Cplex.setParam(IloCplex::StrongItLim, 1);

  Cplex.setParam(IloCplex::WorkDir, ".");
  Cplex.setParam(IloCplex::NodeFileInd, 2);
  //TreLim

  Cplex.setParam(IloCplex::HeurFreq, -1); // heuristic frequency
  Cplex.setParam(IloCplex::RINSHeur, -1); // do not apply RINS heuristic
  Cplex.setParam(IloCplex::FPHeur, -1); // do not apply feasibility pump heuristic
  Cplex.setParam(IloCplex::LBHeur, 0); // do not apply local branching heuristic
  Cplex.setParam(IloCplex::PreInd, 0); // do not apply presolve
  Cplex.setParam(IloCplex::PreslvNd, -1); // node presolve
  Cplex.setParam(IloCplex::Symmetry, 0); // symmetry breaking
  Cplex.setParam(IloCplex::AggInd, 0); // do not use any aggregator
  Cplex.setParam(IloCplex::BndStrenInd, 0); // no var bound strengthening
  Cplex.setParam(IloCplex::CoeRedInd, 0); // no coefficient reduction
  Cplex.setParam(IloCplex::DepInd, 0); // no dependency checker
  Cplex.setParam(IloCplex::Reduce, CPX_PREREDUCE_NOPRIMALORDUAL); // no reductions
  Cplex.setParam(IloCplex::CutPass, -1); // cutting plane passes at root node
  Cplex.setParam(IloCplex::Cliques, -1);
  Cplex.setParam(IloCplex::Covers, -1);
  Cplex.setParam(IloCplex::DisjCuts, -1);
  Cplex.setParam(IloCplex::FlowCovers, -1);
  Cplex.setParam(IloCplex::FlowPaths, -1);
  Cplex.setParam(IloCplex::FracCuts, -1);
  Cplex.setParam(IloCplex::GUBCovers, -1);
  Cplex.setParam(IloCplex::ImplBd, -1);
  Cplex.setParam(IloCplex::MIRCuts, -1);
  Cplex.setParam(IloCplex::MCFCuts, -1);
  Cplex.setParam(IloCplex::ZeroHalfCuts, -1);

  Cplex.setParam(IloCplex::RandomSeed, 13131313);

}

int TSPPDMS::readData(const char name[]) {
  FILE* inputFile;
  char* pch;
  char str[50];
  int idx;
  Arc* e;

  int u, v, w;
  int dummy;
  float xcoord, ycoord;

  if ((inputFile = fopen(name, "r")) == NULL)
    return 1;

  if (fgets(str, 50, inputFile) == NULL)
    return 1;

  //reads the number of requsitions n
  sscanf(str, "%d ", &n);
  V = 2 * n + 2; //the number of nodes 0...2n+1
  A = (V * (V - 1)) / 2; //number of edges (complete graph)

  memset(is_pickup, false, sizeof (is_pickup));


  //reads the points (one for each line, in the next nVertices lines)
  while (nodesCoord.size() != V) {
    if (fgets(str, 50, inputFile) != NULL) {
      sscanf(str, "%d %f %f", &dummy, &xcoord, &ycoord);
      nodesCoord.push_back(make_pair(xcoord, ycoord));
    } else
      return 1;
  }

  //read the pickup and delivery assigments
  for (int i = 0; i < n; i++) {
    if (fgets(str, 50, inputFile) != NULL) {
      sscanf(str, "%d %d %d", &u, &w, &v);
      is_pickup[u] = true; //node u is a pickup
      //      is_pickup[v] = false; //node v is a delivery
      sibling[u] = v;
      sibling[v] = u;
      pickups.push_back(u);
      deliveries.push_back(v);
      itemLength[u] = w;
      itemLength[v] = -w;
    } else
      return 1;
  }
  is_pickup[0] = true; //node 0 is the initial depot
  sibling[0] = 0;
  sibling[2 * n + 1] = 2 * n + 1;

  float p1, p2;
  for (u = 0; u < V; u++) {
    for (v = u + 1; v < V; v++) {
      e = new Arc;
      e->u = u;
      e->v = v;
      p1 = nodesCoord[u].first - nodesCoord[v].first;
      p2 = nodesCoord[u].second - nodesCoord[v].second;
      e->c = (int) (sqrt(p1 * p1 + p2 * p2) + 0.5);
      arcs.push_back(e);
      adj[u].push_back(e);
      adj[v].push_back(e);
    }
  }

  //reads the number of available stacks and the capacity
  if (fgets(str, 50, inputFile) != NULL) {
    sscanf(str, "%d %d", &K, &Q);
  } else
    return 1;

  //read the initial upper bound
  if (fgets(str, 50, inputFile) != NULL) {
    sscanf(str, "%d", &initial_ub);
  } else
    return 1;

  filename = string(name);
  //  filename = filename.substr(7); // removing ./inst/

#ifdef debug
  printf("Vertices:%d Edges: %d Requisitions %d\n", V, A, n);
  printf("pickups:\t");
  for (int i = 0; i < pickups.size(); i++)
    printf("%d\t", pickups[i]);
  printf("\n");
  printf("deliveries:\t");
  for (int i = 0; i < deliveries.size(); i++)
    printf("%d\t", deliveries[i]);
  printf("\n");
  printf("Items length:");
  for (int i = 0; i < pickups.size(); i++)
    printf("(%d %d)=(%d %d)\t", pickups[i], sibling[pickups[i]], itemLength[pickups[i]], itemLength[sibling[pickups[i]]]);
  printf("\n");
#endif
  fclose(inputFile);

  return 0;
}

void TSPPDMS::create_graphviz_image(int n, int solution_type) {
  if (solution_type == INTEGER_SOLUTION) {
    char command[255];
    int result;
    sprintf(command, "neato -Tgif ./graphs/graph.txt > ./graphs/graph%d.gif", n);
    create_graphviz_solution();
    result = system(command);
    result = system("rm ./graphs/graph.txt");
  } else if (solution_type == LINEAR_SOLUTION) {
    char command[255];
    int result;
    sprintf(command, "neato -Tgif ./graphs/graph.txt > ./graphs/graph%d.gif", n);
    create_graphviz_linear_solution();
    result = system(command);
    result = system("rm ./graphs/graph.txt");
  }
}

/*graph G
{
        n1 [pos = "0,1!", height = 0.3, label ="", color = "red"];

        n1 -- n2 [dir = forward, color = blue];

}*/

void TSPPDMS::create_graphviz_linear_solution() {
  const char* colors[] = {"aquamarine", "blue", "bisque1", "burlywood", "cadetblue", "chartreuse",
    "coral", "cornflowerblue", "crimson", "darkorange", "darksalmon", "darkviolet", "gold", "goldenrod",
    "khaki", "greenyellow", "lightgrey", "palegoldenrod", "palegreen", "orangered", "tomato", "pink", "plum"};
  ofstream myfile("./graphs/graph.txt");
  double scale;
  if (nodesCoord[0].first > 10000) {
    scale = 1000;
  } else if (nodesCoord[0].first > 1000) {
    scale = 100;
  } else {
    scale = 10;
  }
  if (myfile.is_open()) {
    myfile << "graph G" << "\n" << "{\noverlap=false;\n";
    myfile << "n" << 0 << "[pos = \"" << nodesCoord[0].first / scale << "," << nodesCoord[0].second / scale
            << "\", height = 0.3, width = 0.3, label =\"" << 0 << "\", fixedsize=\"true\"];\n";
    for (int i = 0; i < pickups.size(); i++) {
      myfile << "n" << pickups[i] << "[pos = \"" << nodesCoord[pickups[i]].first / scale << "," << nodesCoord[pickups[i]].second / scale
              << "\", height = 0.3, width = 0.3, label =\"" << pickups[i] << "\", style=\"filled\", fillcolor = \"" << colors[i] << "\";fixedsize=\"true\"];\n";
      myfile << "n" << sibling[pickups[i]] << "[pos = \"" << nodesCoord[sibling[pickups[i]]].first / scale << "," << nodesCoord[sibling[pickups[i]]].second / scale
              << "\", height = 0.3, width = 0.3, shape=\"doublecircle\", label =\"" << sibling[pickups[i]] << "\", style=\"filled\", fillcolor = \"" << colors[i] << "\";fixedsize=\"true\"];\n";
    }
    myfile << "n" << 2 * n + 1 << "[pos = \"" << nodesCoord[2 * n + 1].first / scale << "," << nodesCoord[2 * n + 1].second / scale
            << "\", height = 0.3, width = 0.3, shape=\"doublecircle\", label =\"" << 2 * n + 1 << "\", fixedsize=\"true\"];\n";
    for (int i = 0; i < pickups.size(); i++) {
      if (_x[mapVarsX[100 * 0 + pickups[i]]] > 0) {
        myfile << "n" << 0 << "-- n" << pickups[i] << " [dir = \"forward\", color = black];\n";
        for (int k = 0; k < 1; k++) {
          if ((_y[mapVarsY[10 * pickups[i] + k]] > 0)) {
            myfile << "n" << 0 << "-- n" << pickups[i] << " [dir = \"forward\", label=\"k" << k << " = " << _y[mapVarsY[10 * pickups[i] + k]] << "\", color = blue];\n";
          }
        }
      }
    }
    for (int i = 1; i < 2 * n + 1; i++) {
      for (int j = 1; j <= 2 * n + 1; j++) {
        if (i != j && !(!is_pickup[i] && sibling[i] == j) && !(is_pickup[i] && j == 2 * n + 1) && (_x[mapVarsX[100 * i + j]] > 0.01)) {
          myfile << "n" << i << "-- n" << j << " [dir = \"forward\", label=\"" << _x[mapVarsX[100 * i + j]] << "\", color = black];\n";
        }
      }
    }

    for (int j = 0; j < pickups.size(); j++) {
      for (int i = 1; i < 2 * n + 1; i++) {
        if (i != pickups[j] && !(!is_pickup[i] && sibling[i] == pickups[j])) {
          if (_x[mapVarsX[100 * i + pickups[j]]] > 0) {
            for (int k = 0; k < K; k++) {
              if (_y[mapVarsY[10 * pickups[j] + k]] > 0.01) {
                myfile << "n" << i << "-- n" << pickups[j] << " [dir = \"forward\", label=\"k" << k << " = " << _y[mapVarsY[10 * pickups[j] + k]] << "\", color = blue];\n";
              }
            }
            i = 2 * n + 1;
          }
        }
      }
    }

    myfile << "}";
  }
  myfile.close();
}

void TSPPDMS::create_graphviz_solution() {
  const char* colors[] = {"aquamarine", "blue", "bisque1", "burlywood", "cadetblue", "chartreuse",
    "coral", "cornflowerblue", "crimson", "darkorange", "darksalmon", "darkviolet", "gold", "goldenrod",
    "khaki", "greenyellow", "lightgrey", "palegoldenrod", "palegreen", "orangered", "tomato", "pink", "plum"};
  ofstream myfile("./graphs/graph.txt");
  double scale;
  if (nodesCoord[0].first > 10000) {
    scale = 1000;
  } else if (nodesCoord[0].first > 1000) {
    scale = 100;
  } else {
    scale = 10;
  }
  if (myfile.is_open()) {
    myfile << "graph G" << "\n" << "{\noverlap=false;\n";
    myfile << "n" << 0 << "[pos = \"" << nodesCoord[0].first / scale << "," << nodesCoord[0].second / scale
            << "\", height = 0.3, width = 0.3, label =\"" << 0 << "\", fixedsize=\"true\"];\n";
    for (int i = 0; i < pickups.size(); i++) {
      myfile << "n" << pickups[i] << "[pos = \"" << nodesCoord[pickups[i]].first / scale << "," << nodesCoord[pickups[i]].second / scale
              << "\", height = 0.3, width = 0.3, label =\"" << pickups[i] << "\", style=\"filled\", fillcolor = \"" << colors[i] << "\";fixedsize=\"true\"];\n";
      myfile << "n" << sibling[pickups[i]] << "[pos = \"" << nodesCoord[sibling[pickups[i]]].first / scale << "," << nodesCoord[sibling[pickups[i]]].second / scale
              << "\", height = 0.3, width = 0.3, shape=\"doublecircle\", label =\"" << sibling[pickups[i]] << "\", style=\"filled\", fillcolor = \"" << colors[i] << "\";fixedsize=\"true\"];\n";
    }
    myfile << "n" << 2 * n + 1 << "[pos = \"" << nodesCoord[2 * n + 1].first / scale << "," << nodesCoord[2 * n + 1].second / scale
            << "\", height = 0.3, width = 0.3, shape=\"doublecircle\", label =\"" << 2 * n + 1 << "\", fixedsize=\"true\"];\n";
    for (int i = 0; i < pickups.size(); i++) {
      for (int k = 0; k < 1; k++) {
        if ((_y[mapVarsY[1000 * 0 + 10 * pickups[i] + k]] > 1 - 0.01) && (_x[mapVarsX[100 * 0 + pickups[i]]] > 1 - 0.01)) {
          myfile << "n" << 0 << "-- n" << pickups[i] << " [dir = \"forward\", label=\"" << k << "\", color = blue];\n";
          myfile << "n" << 0 << "-- n" << pickups[i] << " [dir = \"forward\", color = black];\n";
        }
      }
    }
    for (int i = 1; i < 2 * n + 1; i++) {
      for (int j = 1; j <= 2 * n + 1; j++) {
        if (i != j && !(!is_pickup[i] && sibling[i] == j) && !(is_pickup[i] && j == 2 * n + 1) && (_x[mapVarsX[100 * i + j]] > 1 - 0.01)) {
          myfile << "n" << i << "-- n" << j << " [dir = \"forward\", color = black];\n";
        }
      }
    }
    for (int j = 0; j < pickups.size(); j++) {
      for (int i = 1; i < 2 * n + 1; i++) {
        if (i != pickups[j] && !(!is_pickup[i] && sibling[i] == pickups[j])) {
          if (_x[mapVarsX[100 * i + pickups[j]]] > 1 - EPSILON) {
            for (int k = 0; k < K; k++) {
              if (_y[mapVarsY[10 * pickups[j] + k]] > 1 - EPSILON) {
                myfile << "n" << i << "-- n" << pickups[j] << " [dir = \"forward\", label=\"k" << k << " = " << _y[mapVarsY[10 * pickups[j] + k]] << "\", color = blue];\n";
              }
            }
            i = 2 * n + 1;
          }
        }
      }
    }
    myfile << "}";
  }
  myfile.close();
}

bool TSPPDMS::checkSolution() {
  bool flag = false;
  int u = 0;
  int v;
  int di;
  vector<int> tour;
  int visited[50];
  stack<int> stk[MAXSTKS];
  int maxCap = 0;
  int used_stack[MAXNODES];
  memset(used_stack, 0, sizeof (used_stack));

  tour.push_back(u);
  while (!flag) {
    for (int j = 1; j <= 2 * n + 1; j++) {
      if (u != j && !(!is_pickup[u] && sibling[u] == j) && u != 2 * n + 1
              && !(is_pickup[u] && j == 2 * n + 1) && !(u == 0 && !is_pickup[j])) {
        if (_x[ mapVarsX[100 * u + j] ] > 1 - EPSILON) {
          if (is_pickup[j]) {
            for (int k = 0; k < K; k++) {
              if (_y[ mapVarsY[10 * j + k] ] > 1 - EPSILON) {
                used_stack[j] = k;
                used_stack[sibling[j]] = k;
                k = K;
              }
            }
          }
          u = j;
          j = 2 * n + 2;
          if (u == 2 * n + 1)
            flag = true;
          tour.push_back(u);
        }
      }
    }
  }
  //  printf("%d\n", 2 * n + 1);

  for (int i = 0; i <= 2 * n + 1; i++)
    visited[i] = false;
  flag = true;
  for (int i = 1; i < tour.size() - 1 && flag; i++) {//check if each pickup node is visited before its corresponding delivery node in the solution
    u = tour[i];
    v = sibling[u];
    visited[u] = true;
    if (is_pickup[u] && visited[v]) {
      printf("infeasible solution...\n %d(delivery) visited before %d(pickup)\n", v, u);
      flag = false;
    }
  }

  int k;
  if (flag) {//check if lifo loading policy is respected in the solution
    flag = true;
    u = 0;
    for (int i = 1; i < tour.size() - 1; i++) {
      v = tour[i];
      k = used_stack[v];
      if (is_pickup[v]) {//if its a pickup, load it in the corresponding stack
        stk[k].push(sibling[v]);
        if (maxCap < stk[k].size())
          maxCap = stk[k].size();
      }
      if (!is_pickup[v]) {//if its a delivery, checks if its on the top of the stack
        di = stk[k].top();
        if (di != v) {//shows all crossing itens loaded on a stack...
          printf("infeasible solution... delivery %d crosses delivery %d on stack %d\n", di, v, k);
          flag = false;
        } else
          stk[k].pop();
      }
      u = v;
    }
  }

  int start[MAXSTKS];
  int end[MAXSTKS];
  int capacities[MAXSTKS];
  int old_capacity;
  memset(capacities, 0, sizeof (capacities));
  memset(visited, false, sizeof (visited));
  for (int i = 1; i < tour.size(); i++) {
    v = tour[i];
    k = used_stack[v];
    if (is_pickup[v]) {
      old_capacity = capacities[k];
      capacities[k] += itemLength[v];
      visited[v] = true;
      if (old_capacity == 0) //before loading v, stack k was empty...
      {
        start[k] = i - 1;
      } else if (capacities[k] > Q) {
        end[k] = i;
        printf("Stack capacity is %d but the solution uses %d!\n", Q, capacities[k]);
        flag = false;
      }
    }
    if (!is_pickup[v] && visited[sibling[v]]) {
      capacities[k] -= -itemLength[v];
    }
    u = v;
  }

  return flag;
}
