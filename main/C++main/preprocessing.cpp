/**
preprocessing.py
Before solving an PCI instance using our models, two reduction rules are applied
to the instance, which may generate separated graph components, each component a
new instance.
**/

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#define REMOVED -1
#define INFECTED -2
#define DEGREE_ONE_REMOVAL -3
#define FREE 0
#define MIN_W 1
#define MAX_W 100

template<class T>
void print_vector2(vector<T> &v, string message) {
  cout << message;
  for (int i = 0; i < v.size(); i ++)
    cout << v[i] << "  ";
  cout << endl;
}

template<class T>
void print_matrix2(vector<vector<T> > &m, string message) {
  cout << message;
  for (int i = 0; i < m. size(); i ++)
    print_vector(m[i], "");
}

// dimacs are the coloring problem instances, and they have one extra argument
// where it defines which way to create the trashold vector f
void read_file_dimacs(string filename, string option, vector<vector<int> > 
  &component) {
  cout << "dimacs" << endl;
  int N = 0, M = 0;
  vector<vector<bool> > adjacency_matrix;
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
	adjacency_matrix = vector<vector<bool> >(N, vector<bool>(N, false));
	int u, v;
	for (int i = 0; i < M; i ++) {
		inputfile >> line >> u >> v;
		adjacency_matrix[u - 1][v - 1] = true;
		adjacency_matrix[v - 1][u - 1] = true;
	}
  inputfile.close();
  
  vector<vector<int> > adjacency_list(N);
  for (int i = 0 ; i < N; i ++)
  	for (int j = i + 1; j < N; j ++) {
  		adjacency_list[i].push_back(j);
  		adjacency_list[j].push_back(i);
  	}
  	
  vector<int> f(N), w(N);
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

  adjacency_list.push_back(f);
  adjacency_list.push_back(w);
  component = adjacency_list;
}

// instances created by us are in this format
void read_file_formated(string filename, vector<vector<int> > & component) {
  fstream inputfile(filename);
  int N, M;
  char ch;
  inputfile >> ch;
  if (! inputfile.is_open()) {
    cout << "Unable to open file" << endl; 
    exit(0);
  } 
  else if (ch == 'c') {
    cout << "Dimacs file with too few arguments" << endl; 
    exit(0);
  }
  N = (int) ch;
  inputfile >> M;

  vector<int> f(N), w(N);
  for (int i = 0; i < N; i ++)
    inputfile >> f[i] >> w[i];
  vector<vector<bool> > adjacency_matrix(N, vector<bool>(N, false));
  int u, v;
  for (int i = 0; i < M; i ++) {
    inputfile >> u >> v;

    adjacency_matrix[u][v] = true;
    adjacency_matrix[v][u] = true;
  }
  inputfile.close();
  vector<vector<int> > adjacency_list(N);
  for (int i = 0; i < N; i ++) {
    for (int j = i + 1; j < N; j ++) {
      if (adjacency_matrix[i][j]) {
        adjacency_list[i].push_back(j);
        adjacency_list[j].push_back(i);
      }
    }
  }
  adjacency_list.push_back(f);
  adjacency_list.push_back(w);
  component = adjacency_list;
}


// This reduction infects vertices that have to be infected in the beginning of
//the process: if N_G(v) < f[v], v is in the solution
bool reduction_one(vector<vector<int> > &adjacency_list, vector<int> &f,
  vector<int> &num_neighbors, vector<bool> &removed_vertices,
  vector<int> &type) {
  
  bool infected_someone = false, new_round_needed = true;
  while (new_round_needed) {
    new_round_needed = false;
    for (int v = 0; v < f.size(); v ++) {
      if (! removed_vertices[v] && (f[v] <= 0 || num_neighbors[v] < f[v])) {
        new_round_needed = true;
        infected_someone = true;
        removed_vertices[v] = true;
        if (f[v] <= 0)
          type[v] = REMOVED;
        else
          type[v] = INFECTED;

        for (vector<int>::iterator u = (adjacency_list[v]).begin();
          u != (adjacency_list[v]).end(); ++ u) {
          num_neighbors[* u] -= 1;
          f[* u] -= 1;
        }
      }
    }

  }
  return infected_someone;
}


bool reduction_two (vector<int> &f, vector<int> num_neighbors, 
  vector<bool> &removed_vertices, vector<int> &type) {
  bool removed_vertex = false;
  for (int v = 0; v < f.size(); v ++) {
    if (! removed_vertices[v] && num_neighbors[v] == 1 && f[v] == 1) {
      removed_vertex = true;
      removed_vertices[v] = true;
      type[v] = DEGREE_ONE_REMOVAL;
    } 
  }
  return removed_vertex;
}

vector<vector<int> > reduced_instance (vector<vector<int> > adjacency_list, 
  vector<int> f, vector<int> w, vector<bool> &removed_vertices) {
  int N = f.size();
  vector<vector<bool> > adjacency_matrix(N, vector<bool>(N, false));
  for (int v = 0; v < N; v ++) {
    for (int j = 0; j < adjacency_list[v].size(); j ++) {
      int u = adjacency_list[v][j];
      adjacency_matrix[u][v] = true;
      adjacency_matrix[v][u] = true;
    }
  }
  adjacency_list = vector<vector<int> >();
  vector<int> new_indices;
  int i = 0;
  for (int v = 0; v < N; v ++) {
    if (! removed_vertices[v]) {
      assert(f[v] > 0);
      f[i] = f[v];
      w[i] = w[v];
      new_indices.push_back(v);
      adjacency_list.push_back(vector<int>());
      int j = 0;
      for (int u = 0; u < N; u ++) {
        if (! removed_vertices[u]) {
          if (adjacency_matrix[u][v])
            adjacency_list[i].push_back(j);
          j ++;
        }
      }
      i ++;
    }
  }
  f.resize(i);
  w.resize(i);
  adjacency_list.push_back(f);
  adjacency_list.push_back(w);
  return adjacency_list;
}

deque<vector<vector<int> > >
separate_in_connected_instances(vector<vector<int> > &adjacency_list, 
  vector<int> &f, vector<int> &w) {
  int N = f.size();
  vector<bool> visited(N, false);
  deque<vector<vector<int> > > components;
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
      reduced_instance(adjacency_list, f, w, not_in_component));
  } 
  return components;
}


deque<vector<vector<int> > > data_preprocessing(int argc, char* argv[]) {
  // reads graphs adjacency list, f and w from file
  vector<vector<int> > adjacency_list;
  if (argc == 4) {
    //in case the file is in the dimacs format
    read_file_dimacs(argv[2], argv[3], adjacency_list);
  }
  else if (argc == 3) {
    // file formated by us
    read_file_formated(argv[2], adjacency_list);
  }
  else {
    cout << "PCI_problem_cplex <model> <datafile> <?option if in dimacs format>"
      << endl;
    exit(0);
  }
  vector<int> f, w;
  f = adjacency_list[adjacency_list.size() - 2];
  w = adjacency_list[adjacency_list.size() - 1];
  adjacency_list.pop_back();
  adjacency_list.pop_back();
  int N = f.size();
  for (int v = 0; v < N; v ++) {
    cout << ", " << v;
    if (adjacency_list[v].size() < f[v]) cout << " WRONG ";
  }
  cout << endl;
  vector<bool> removed_vertices(N, false);
  vector<int> num_neighbors(N, 0), type_vertex(N, FREE);
  for (int i = 0; i < N; i ++)
    num_neighbors[i] = adjacency_list[i].size();

  bool vertices_were_removed = true;
  while (vertices_were_removed) {

    vertices_were_removed = reduction_one(adjacency_list, f, num_neighbors,
      removed_vertices, type_vertex);

    vertices_were_removed |= reduction_two(f, num_neighbors,
      removed_vertices, type_vertex);
  }

  for (int v = 0; v < N; v ++) 
    assert(removed_vertices[v] == (type_vertex[v] < 0));
  adjacency_list = reduced_instance(adjacency_list, f, w, removed_vertices);
  f = adjacency_list[adjacency_list.size() - 2];
  w = adjacency_list[adjacency_list.size() - 1];
  adjacency_list.pop_back();
  adjacency_list.pop_back();
  for (int v = 0; v < N; v ++) {
    cout << ", " << v;
    if (adjacency_list[v].size() < f[v]) cout << " WRONG ";
  }
  cout << endl;
  return separate_in_connected_instances(adjacency_list, f, w);
}