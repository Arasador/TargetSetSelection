/*#include <assert.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <deque>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdarg>
#include <random>
#include <utility>      // std::pair, std::make_pair

using namespace std;
/*/
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
#define INCLUDE

enum removals { REMOVED = -1, INFECTED = -2, DEGREE_ONE_REMOVAL = -3, FREE = 0};

enum model {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
  S_SMALLER_H1, S_SMALLER_H2, S_SMALLER_NEW};

// list of all models we have
  model models[] = {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
  S_SMALLER_H1, S_SMALLER_H2, S_SMALLER_NEW};


#define TIMELIMIT 180000
//#define PRINT_LOG
//#define FILE_S_CUTTER_INFO
#define ROOT_RELAX
#define UPPERBOUND_CUT
#define USERCUT
#define EXCLUDE_SMALL_CLUSTERS
#define MIN_SIZE_CLUSTER 3
//#define CONSIDER_WEIGHTS
#define CONSIDER_F_RANDOM