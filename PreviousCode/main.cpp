#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
//#include <omp.h>
#include <ilcplex/ilocplex.h>
#include "TSPPDMS.h"

using namespace std;

int main(int argc, char** argv) {
  int alg;
  const char* fileName = NULL;
  const char* svgName = NULL;
  const char* lpName = NULL;
  const char* parName = NULL;
  char param;
  int numThreads = 1;
  int vlEmph = 0;
  int vlDisp = 0;
  double vlGap = 0.0;
  int memory = 0;

  bool cutFound;
  int iter;
  int NumSecs, NumPrecs, NumLifos, NumCaps;

  while ((param = getopt(argc, argv, "f:a:s:m:l:t:g:e:v:h?")) != -1)
    switch (param) {
      case 'f'://[data file]
        fileName = optarg;
        break;
      case 's'://[svg tour file]
        svgName = optarg;
        break;
      case 'a'://[simplex algorithm]
        alg = atoi(optarg);
        break;
      case 'l'://[lp file model]
        lpName = optarg;
        break;
      case 'm'://[CPLEX working memory]
        memory = atoi(optarg);
        break;

      case 't'://numero de threads
        sscanf(optarg, "%d", &numThreads);
        break;
      case 'g'://EPgap
        sscanf(optarg, "%lf", &vlGap);
        break;
      case 'e'://MIPEmphasis
        sscanf(optarg, "%d", &vlEmph);
        break;
      case 'v'://CPLEX verbose
        sscanf(optarg, "%d", &vlDisp);
        break;
      case 'h': case '?'://help
        printf("usage: ./pdtspms -f [input_file] ...\n");
        exit(0);
        break;
      default:
        printf("invalid argument...!\n");
    }

  if (fileName == NULL) {
    cout << "Sem arquivo para ler!\n";
    return 1;
  }

  IloEnv env_cp;
  TSPPDMS tsppdms(env_cp);
  try {
    if (tsppdms.readData(fileName) == 1) {
      cout << "error while reading the input file...!\n";
      return 1;
    }
    tsppdms.setCplexSettings(vlDisp, vlEmph, alg, numThreads, vlGap, memory);

    tsppdms.setModelProblem();
    tsppdms.startAlg();

    tsppdms.relaxIntVars();

    cutFound = true;
    iter = 0;
    tsppdms.startAlg();
    while (cutFound) {
      tsppdms.solveProblem();
      cutFound = tsppdms.findCuts(NumSecs, NumPrecs, NumLifos,NumCaps);
      if (cutFound) {
        //printf("%d | %.2f | %d | %d | %d | %d |\n", ++iter, tsppdms.getSolution(), NumSecs, NumPrecs, NumLifos, NumCaps);
        tsppdms.addCuts();
      }
    }
    
    tsppdms.enforceIntVars();
    tsppdms.solveProblem();
    tsppdms.endAlg();

    if (tsppdms.checkSolution()) {
//      printf("solution checked!\n");
    }
  } catch (IloException& ex) {
    tsppdms.Texception();
  }

  return 0;
}