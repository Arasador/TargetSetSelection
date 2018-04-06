//#ifndef SCUTTER_H
//#define SCUTTER_H
#include <assert.h>
#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include <stack>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdarg>
#include <ilcplex/ilocplex.h>
#include <random>
#include <utility>      // std::pair, std::make_pair


using namespace std;

#define S_MODEL 1
#define S_SMALLER 2
#define WS_SMALLER 3
#define DOMINATED 4
#define WDOMINATED 5
#define S_SMALLER_H1 6
#define S_SMALLER_H2 7

#define INFECTED -2

typedef vector<vector<int> > Vector2_int;
typedef vector<vector<bool> > Vector2_bool;
typedef vector<int>::iterator Vector_it;
