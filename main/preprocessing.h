/**
  preprocessing.h
  Before solving an PCI instance using our models, reduction rules are applied
  to the instance, which may generate a new graph and separated graph components,
  each component a new instance.
**/
#ifndef PREPROCESSING
#define PREPROCESSING

#include <string>
#include <deque>
#include <vector>

#include "includes.h"


//enum removals { REMOVED = -1, INFECTED = -2, DEGREE_ONE_REMOVAL = -3, FREE = 0};

using namespace std;

template<class T>
void print_vector2(vector<T>& v, string message);

template<class T>
void print_matrix2(vector<vector<T>>& m, string message);

void read_file_dimacs(string filename, string option, vector<vector<int>>&
  adjacency_list, vector<int>& f, vector<int>& w);

void read_file_formated(string filename, vector<vector<int>>& adjacency_list,
  vector<int>& f, vector<int>& w);

bool reduction_one(const vector<vector<int>>& adjacency_list, vector<int>& f, 
	vector<int>& num_neighbors, vector<bool>& removed_vertices, 
	vector<removals>& type);

bool reduction_two ( vector<vector<int>>& adjacency_list, vector<int>& f,
	vector<int>& num_neighbors, vector<bool>& removed_vertices, 
	vector<removals>& type);

void reduced_instance (vector<vector<int>>& adjacency_list, vector<int>& f, 
  vector<int>& w, vector<bool>& removed_vertices);

vector<vector<int>> reduce_and_return_component(
  vector<vector<int>> adjacency_list, vector<int> f, vector<int> w,
  vector<bool>& not_in_component);

deque<vector<vector<int>>>
separate_in_connected_instances(vector<vector<int>> &adjacency_list, 
  vector<int> &f, vector<int> &w);

deque<vector<vector<int>>> data_preprocessing(int argc, char* argv[]);

#endif  