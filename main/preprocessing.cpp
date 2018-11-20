/**
  preprocessing.cpp
**/
 
#include "preprocessing.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <assert.h>

//#include "includes.h"
/*#define REMOVED -1
#define INFECTED -2
#define DEGREE_ONE_REMOVAL -3
#define FREE 0 //*/
#define MIN_W 1
#define MAX_W 100

//using namespace std;


template<class T>
void print_vector2(vector<T>& v, string message) {
  cout << message;
  for (int i = 0; i < v.size(); i ++)
    cout << v[i] << "  ";
  cout << endl;
}

template<class T>
void print_matrix2(vector<vector<T>>& m, string message) {
  cout << message;
  for (int i = 0; i < m. size(); i ++)
    print_vector(m[i], "");
}

// dimacs are the coloring problem instances, and they have one extra argument
// where it defines which way to create the trashold vector f
void read_file_dimacs(string filename, string option, vector<vector<int>>&
  adjacency_list, vector<int>& f, vector<int>& w) {
  cout << "dimacs" << endl;
  int N = 0, M = 0;
  vector<vector<bool>> adjacency_matrix;
  ifstream inputfile (filename);
  if (! inputfile.is_open()) {
    cout << "Unable to open file"; 
    abort();
  }
	string id, line;
	inputfile >> id;
	while (id.compare("c") == 0) {
	  getline (inputfile,line);
	  inputfile >> id;
	}
	inputfile >> line >> N >> M;
	adjacency_matrix = vector<vector<bool>>(N, vector<bool>(N, false));
	int u, v;
	for (int i = 0; i < M; i ++) {
		inputfile >> line >> u >> v;
		adjacency_matrix[u - 1][v - 1] = true;
		adjacency_matrix[v - 1][u - 1] = true;
	}
  inputfile.close();
  
  adjacency_list = vector<vector<int>>(N);
  for (int i = 0 ; i < N; i ++)
  	for (int j = i + 1; j < N; j ++) {
  		adjacency_list[i].push_back(j);
  		adjacency_list[j].push_back(i);
  	}
  	
  f = vector<int>(N);
  w = vector<int>(N);
  for (int i = 0; i < N; i ++) {
    w[i] = rand() % (MAX_W - MIN_W) + MIN_W;
    if (option.compare("2") == 0)
      f[i] = 2;
    else if (option.compare("degree") == 0)
      f[i] = max((int) adjacency_list[i].size() / 2, 1);
    else if (option.compare("e2v") == 0)
      f[i] = max(M / (2 * N), 2);
    else if (option.compare("random") == 0)
      f[i] = rand() % adjacency_list[i].size() + 1;
  }
}

// instances created by us are in this format
void read_file_formated(string filename, vector<vector<int>>& adjacency_list,
  vector<int>& f, vector<int>& w) {
  fstream inputfile(filename);
  int N, M;
  string first_char_in_file;
  inputfile >> first_char_in_file;
  if (! inputfile.is_open()) {
    cout << "Unable to open file" << endl; 
    exit(0);
  } 
  else if (first_char_in_file.compare("c") == 0) {
    cout << "Dimacs file with too few arguments" << endl; 
    exit(0);
  }
  N = stoi(first_char_in_file);
  inputfile >> M;

  f = vector<int>(N);
  w = vector<int>(N);
  for (int i = 0; i < N; i ++)
    inputfile >> f[i] >> w[i];
  vector<vector<bool>> adjacency_matrix(N, vector<bool>(N, false));
  int u, v;
  for (int i = 0; i < M; i ++) {
    inputfile >> u >> v;
    adjacency_matrix[u][v] = true;
    adjacency_matrix[v][u] = true;
  }
  inputfile.close();
  adjacency_list = vector<vector<int>>(N);
  for (int i = 0; i < N; i ++) {
    for (int j = i + 1; j < N; j ++) {
      if (adjacency_matrix[i][j]) {
        adjacency_list[i].push_back(j);
        adjacency_list[j].push_back(i);
        adjacency_matrix[i][j] = false;
        adjacency_matrix[j][i] = false;
      }
    }
  }
  for (int i = 0; i < N; i ++) {
    for (int j = 0; j < N; j ++) {
      assert(! adjacency_matrix[i][j]);
    }
  }
}


// This reduction infects vertices that have to be infected in the beginning of
//the process: if N_G(v) < f[v], v is in the solution
bool reduction_one(const vector<vector<int>>& adjacency_list, 
  vector<int>& f, vector<int>& num_neighbors, 
  vector<bool>& removed_vertices, vector<removals>& type) {
  
  bool infected_someone = false, new_round_needed = true;
  // while removed someone, has to check if someone was affected by that
  while (new_round_needed) {
    new_round_needed = false;
    // tests if any vertex needs to be infected initially, or is infected 
    // automatically
    for (int v = 0; v < f.size(); v ++) {
      if (! removed_vertices[v] && (f[v] <= 0 || num_neighbors[v] < f[v])) {
        new_round_needed = true;
        infected_someone = true;
        removed_vertices[v] = true;
        if (f[v] <= 0)
          type[v] = REMOVED;
        else
          type[v] = INFECTED;

        for (auto u: adjacency_list[v]) {
          num_neighbors[u] --;
          f[u] -= 1;
        }
      }
    }

  }
  return infected_someone;
}

// If num neighbors is 1 and f[v] == 1, remove vertex, does not affect
// rest of the simulation
bool reduction_two ( vector<vector<int>>& adjacency_list, 
  vector<int>& f, vector<int>& num_neighbors, 
  vector<bool>& removed_vertices, vector<removals>& type) {
  bool removed_vertex = false, new_round_needed = true;
  while (new_round_needed) {
    // needs new round if removed someone (can affect other vertices)
    new_round_needed = false;
    for (int v = 0; v < f.size(); v ++) {
      if (! removed_vertices[v] && num_neighbors[v] == 1 && f[v] == 1) {
        removed_vertex = true;
        removed_vertices[v] = true;
        new_round_needed = true;
        type[v] = DEGREE_ONE_REMOVAL;
        // for that one neighbors, should reduce number of neighbors, others
        // are infected or removed, do for all in case
        for (auto u: adjacency_list[v]) {
          num_neighbors[u] --;
        } 
      } 
    }
  }
  return removed_vertex;
}

void reduced_instance (vector<vector<int>>& adjacency_list, vector<int>& f, 
  vector<int>& w, vector<bool>& removed_vertices) {
  int N = f.size();
  vector<int> num_neighbors(N);
  

  vector<vector<bool>> adjacency_matrix(N, vector<bool>(N, false));
  for (int v = 0; v < N; v ++) {
    for (int j = 0; j < adjacency_list[v].size(); j ++) {
      int u = adjacency_list[v][j];
      adjacency_matrix[u][v] = true;
      adjacency_matrix[v][u] = true;
    }
  }
  vector<int> old_index;
  int i = 0, new_N = 0;
  for (int v = 0; v < N; v ++) {
    if (! removed_vertices[v]) {
      assert(f[v] > 0);
      f[i] = f[v];
      w[i ++] = w[v];
      new_N ++;
      old_index.push_back(v);
    }
  }
  f.resize(new_N);
  w.resize(new_N);
  vector<vector<int>> old_adjacency_list(adjacency_list);
  adjacency_list = vector<vector<int>>(new_N);
  for (int v = 0; v < new_N; v ++) {
    for (int u = v + 1; u < new_N; u ++) {
      // if were neighbors before, are now as well
      if (adjacency_matrix[old_index[u]][old_index[v]]) {
        adjacency_list[v].push_back(u);
        adjacency_list[u].push_back(v);
      }
    }
  }
}

vector<vector<int>> reduce_and_return_component(
  vector<vector<int>> adjacency_list, vector<int> f, vector<int> w,
  vector<bool>& not_in_component) {
  reduced_instance(adjacency_list, f, w, not_in_component);
  for (int v = 0; v < f.size(); v ++) {
    // reduction 1
    assert(adjacency_list[v].size() >= f[v] && adjacency_list[v].size() > 0);
    // reduction 2
    assert(! (adjacency_list[v].size() == 1 && f[v] == 1));
  }
  adjacency_list.push_back(f);
  adjacency_list.push_back(w);
  return adjacency_list;
}

deque<vector<vector<int>>>
separate_in_connected_instances(vector<vector<int>>& adjacency_list, 
  vector<int>& f, vector<int>& w) {
  int N = f.size();
  vector<bool> visited(N, false);
  deque<vector<vector<int>>> components;
  for (int i = 0; i < N; i ++) {
    if (visited[i])
      continue;
    vector<bool> not_in_component(N, true);
    deque<int> queue(1, i);
    while (! queue.empty()) {
      int v = queue[0];
      queue.pop_front();
      visited[v] = true;
      not_in_component[v] = false;
      for (vector<int>::iterator u = (adjacency_list[v]).begin(); u != 
        (adjacency_list[v]).end(); ++ u) {
        if (! visited[* u]) {
          queue.push_back(* u);
          visited[* u] = true;
        }
      }
    }
    

    components.push_back(
      reduce_and_return_component(adjacency_list, f, w, not_in_component)
    );
      
  } 
  return components;
}

deque<vector<vector<int>>> data_preprocessing(int argc, char* argv[]) {
  // reads graphs adjacency list, f and w from file
  vector<vector<int>> adjacency_list;
  vector<int> f, w;
  if (argc == 4) {
    //in case the file is in the dimacs format
    read_file_dimacs(argv[2], argv[3], adjacency_list, f, w);
  }
  else if (argc == 3) {
    // file formated by us
    read_file_formated(argv[2], adjacency_list, f, w);

  }
  else {
    cout << "PCI_problem_cplex <model> <datafile> <?option if in dimacs format>"
      << endl;
    exit(0);
  }
  int N = f.size();
  vector<bool> removed_vertices(N, false);
  vector<int> num_neighbors(N, 0);
  vector<removals> type_vertex(N, FREE);
  for (int i = 0; i < N; i ++)
    num_neighbors[i] = adjacency_list[i].size();
  // while is removing vertices from instance, continue to try to reduce
  // more the problem
  bool vertices_were_removed = true;
  while (vertices_were_removed) {
    // f[v] > N[v] => vertex need to be in the initial infection
    vertices_were_removed = reduction_one(adjacency_list, f, num_neighbors,
      removed_vertices, type_vertex);
    // if  N[v] == f[v] == 1, this vertex will be infected automatically
    vertices_were_removed |= reduction_two(adjacency_list, f, num_neighbors,
      removed_vertices, type_vertex);
  }
  for (int v = 0; v < N; v ++) {
    assert(removed_vertices[v] == (type_vertex[v] < 0));
    assert(num_neighbors[v] >= f[v] || removed_vertices[v]);
  }
  // take out infected vertices to have a new instance
  reduced_instance(adjacency_list, f, w, removed_vertices);
  for (int v = 0; v < f.size(); v ++) {
    // reduction 1
    assert(adjacency_list[v].size() >= f[v] && adjacency_list[v].size() > 0);
    // reduction 2
    assert(! (adjacency_list[v].size() == 1 && f[v] == 1));
  }
  return separate_in_connected_instances(adjacency_list, f, w);
}