#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <deque>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdarg>
#include <ilcplex/ilocplex.h>
#include <random>
#include <utility>      // std::pair, std::make_pair

using namespace std;

#define NUM_MODELS 8
/*#define S_MODEL 1
#define S_SMALLER 2
#define WS_SMALLER 3
#define DOMINATED 4
#define WDOMINATED 5
#define S_SMALLER_H1 6
#define S_SMALLER_H2 7
// #define S_SMALLER_NEW 8 */
//#define INFECTED -2

enum removals { REMOVED = -1, INFECTED = -2, DEGREE_ONE_REMOVAL = -3, FREE = 0};

enum model {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
  S_SMALLER_H1, S_SMALLER_H2, S_SMALLER_NEW};

// list of all models we have
  model models[] = {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
  S_SMALLER_H1, S_SMALLER_H2, S_SMALLER_NEW};

typedef vector<vector<int> > vectorint2;
typedef vector<vector<bool> > vectorbool2;
//typedef vector<int>::iterator vector_it;

#define TIMELIMIT 30000
//#define PRINT_LOG
//#define FILE_S_CUTTER_INFO
//#define USERCUT
#define EXCLUDE_SMALL_CLUSTERS
#define MIN_SIZE_CLUSTER 3
