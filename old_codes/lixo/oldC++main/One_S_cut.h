#ifndef SCUTTER_H
#define SCUTTER_H

#include"S_cutter.h"
#include <iostream>
#include <vector>
#include <map>
#include <stack>
#include <algorithm>
#include <cstdarg>
#include<string.h>
#include <ilcplex/ilocplex.h>

class One_S_cut: public S_cutter {
public:
  One_S_cut(Vector2_int, vector<int>);
  ~One_S_cut();
}

One_S_cut::One_S_cut(Vector2_int _adjacency_list, vector<int> _f) {
  adjacency_list = _adjacency_list;
	f = _f;
	N = f.size();
}
One_S_cut::~One_S_cut() {}

//------ FIRST SIMPLE S MODEL

// gets the not infected vertices, uses that as s and creates constraints
bool One_S_cut::finds_s_model_constraints(vector<bool> infected_vertices,
  string option) {
    vector<int> new_f = infect_graph(infected_vertices);
    return finds_S_constraints(infected_vertices, new_f);
}
