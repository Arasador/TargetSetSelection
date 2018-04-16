
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
#define S_MODEL 1
#define S_SMALLER 2
#define WS_SMALLER 3
#define DOMINATED 4
#define WDOMINATED 5
#define S_SMALLER_H1 6
#define S_SMALLER_H2 7
#define S_SMALLER_NEW 8
#define INFECTED -2

int models[] = {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
  S_SMALLER_H1, S_SMALLER_H2, S_SMALLER_NEW};

typedef vector<vector<int> > Vector2_int;
typedef vector<vector<bool> > Vector2_bool;
typedef vector<int>::iterator Vector_it;
