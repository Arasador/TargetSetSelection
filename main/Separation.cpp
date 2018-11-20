#include "Separation.h"
#include "S_cutter.h"
#include "ClusterHCS.cpp"


Separation::Separation(vector<vector<int>> _adjacency_list,
	vector<int> _f, vector<int> _w) {
	adjacency_list = _adjacency_list;
	f = _f;
	w = _w;
	N = _f.size();
	clusters = HCS(adjacency_list);
}

Separation::~Separation() {}

int Separation::neighbors_outside_s(const vector<int> &S, int v) {
  int num_neighbors_outside = 0;
  vector<bool> selected_in_S(N, false);
  // creates vector with selected s values seted to true
  for (int i = 0 ; i < S.size(); i ++)
    selected_in_S[S[i]] = true;
  for (vector<int>::iterator u = (adjacency_list[v]).begin();
    u != (adjacency_list[v]).end(); ++ u) {
    if (! selected_in_S[* u])
      num_neighbors_outside ++;
  }
  return num_neighbors_outside;
}

int Separation::S_induced_inequality(const vector<int> &S) {
  //assert(! S.empty());
  int S_lower_bound = f[S[0]] - neighbors_outside_s(S, S[0]);
  for (int v = 1; v < S.size(); v ++) {
    int current_value = f[S[v]] - neighbors_outside_s(S, S[v]);
    if (current_value < S_lower_bound)
      S_lower_bound = current_value;
  }
  return S_lower_bound;
}

int Separation::neighbors_outside_s_bool(const vector<bool> &selected_in_S, int v) {
  int num_neighbors_outside = 0;
  for (auto u : adjacency_list[v]) {
    if (! selected_in_S[u])
      num_neighbors_outside ++;
  }
  return num_neighbors_outside;
}

int Separation::S_induced_inequality_bool(const vector<bool> &S) {
  //assert(! S.empty());
  int S_lower_bound = N;
  for (int v = 0; v < N; v ++) {
		if (! S[v]) continue;
    int current_value = f[v] - neighbors_outside_s_bool(S, v);
    if (current_value < S_lower_bound)
      S_lower_bound = current_value;
  }
  return S_lower_bound;
}

/*template <class T>
int get_argmin(vector<T>& v, vector<bool>& selected) {
	int argmin = 0;
	T min = v[0];
	for (int i = 0; i < N; i ++) {
		if (selected[i] && v[i] < min) {
			min = v[i];
			argmin = i;
		}
	}
	return argmin;
	//return distance(begin(v), min_element(begin(v), end(v)));
}
//*/
// for now, selects vertex with min cost


int Separation::initial_vertex_selection(vector<double>& infected,
	 const vector<bool>& cluster) {

	double min = 2; // all numbers must be non negative and between 0 and 1
	int argmin = -1;
	for (int i = 0; i < N; i ++) {
		if (! cluster.empty() && ! cluster[i]) continue;
		if (min > infected[i]) {
			min = infected[i];
			argmin = i;
		}
	}
	assert(argmin != -1);
	return argmin;
	//auto it = min_element(infected.begin(), infected.end());
	//assert(argmin == it -infected.begin() );
	//return it - infected.begin();
}

// adds neighbors of a given vertex v0 as neighbors of S
void Separation::add_to_neighbors(vector<bool>& neighbors_S, int v0,
	const vector<bool>& cluster) {
	neighbors_S[v0] = false;
	for (auto v: adjacency_list[v0]) {
		if (cluster.empty() || cluster[v]) {
			neighbors_S[v] = true;
		}
	}

}

// here we want to find the v in S's neighborhood that maximizes the right side of the constraint:
// c0(S) >= rhs = min_{v in S}(f(v) - N_{V/S}(v)) (this side needs to be the maximum possible)
// at the same time, you want to minimize c0(S), get the vertex that maximizes rhs - c(v)
int Separation::get_max_S_neighbor(const vector<bool> &neighbors, vector<bool>& S,
	int& rhs, const vector<double>& infected) {
	double maxval = -N;
	int argmax = -1;
	int rhs_violated = S_induced_inequality_bool(S);
	for (int v = 0; v < N; v ++) {
		if (neighbors[v] && ! S[v]) {
			// the whole function S_induced_inequality_bool is not needed
			// we just see if adding v will bring rhs down and use this in the maxval
			int rhs = min(rhs_violated, f[v] - neighbors_outside_s_bool(S, v));
			if (rhs - infected[v] > maxval) {
				maxval = rhs_violated - infected[v];
				argmax = v;
			}
		}
	}
	rhs = maxval;
	return argmax;
}

template<class T>
double get_c0(vector<T>& v) {
	double c0 = 0;
	for (auto i: v) {
		c0 += i;
	}
	return c0;
}

void Separation::initial_heuristic(vector<double>& infected,
	const vector<bool>& cluster) {
	// initializes vectors

	vector<bool> neighbors_S(N, false), S(N, false), best_S(0);
	// for now, takes vertex with least cost
	int v0 = initial_vertex_selection(infected, cluster);
	S[v0] = true;
	double c0 = infected[v0], best_c0 = infected[v0];
	int best_rhs = S_induced_inequality_bool(S);
	best_S = S;
	// adds adjacent vertices to neighbors of S, since S = {v0}
	add_to_neighbors(neighbors_S, v0, cluster);

	while (true) {
		int rhs;
		// gets best vertex : max(rhs - c0(S)) (c0 small rhs big)
		int v = get_max_S_neighbor(neighbors_S, S, rhs, infected);
		if (v == -1) { break; }
		// adds v to S
		S[v] = true;
		neighbors_S[v] = false;
		c0 += infected[v];
		// adds neighbors of v to neighborhood of S

    add_to_neighbors(neighbors_S, v, cluster);
		// in case we have the best c0 till now
		if (best_rhs - best_c0 < rhs - c0) {
			best_rhs = rhs; best_c0 = c0; best_S = S;
		}
		// work on function to take out some vertices.
		//see_if_taking_off_vertex_helps(neighbors_S, )
	}
	if (best_rhs > 0 && best_c0 + 0.0001 < best_rhs) {
		cout << "FOOOOOUND ONE!!!!!!!! bestrhs " <<best_rhs << " > bestc0 " << best_c0 << endl;
		vector<int> int_S;
		for (int i = 0; i < N; i ++) {
			if (S[i]) int_S.push_back(i);
		}
	  constraints_lhs_res.push_back(int_S);
	  constraints_rhs_res.push_back(best_rhs);
	}

}



void Separation::check_if_violates(vector<int> S, vector<int>& rhs,
	vector<vector<int>>& lhs, vector<double>& infected) {
	int rhs_violated = S_induced_inequality(S);
	double c0 = 0;
	for (auto i: S) {
		c0 += infected[i];
	}
	if (c0 < rhs_violated) {
		rhs.push_back(rhs_violated);
		lhs.push_back(S);
	}
}

void Separation::generates_Ss_and_tests(int i, vector<int>& S, int N, vector<int>& rhs,
	vector<vector<int>>& lhs, vector<double>& infected, const vector<bool>& cluster) {
	if (i >= N) {
		if (S.empty()) return;
		check_if_violates(S, rhs, lhs, infected);
		return;
	}
	generates_Ss_and_tests(i + 1, S, N, rhs, lhs, infected);
	// if restricting to clusters, if is not in cluster, returns
	if (cluster.empty() || !cluster[i]) return;
	S.push_back(i);
	generates_Ss_and_tests(i + 1, S, N, rhs, lhs, infected);
	S.pop_back();
}

void Separation::tests_all_Ss(vector<double>& infected) {
	vector<int> S, rhs_violated_s;
	vector<vector<int>> lhs_violated_s;
	generates_Ss_and_tests(0, S, infected.size(), rhs_violated_s, lhs_violated_s,
	infected);
}


void Separation::cluster_find_all(vector<double>& infected, const vector<bool>& cluster) {
	vector<int> S, rhs_violated_s;
	vector<vector<int>> lhs_violated_s;
	generates_Ss_and_tests(0, S, infected.size(), rhs_violated_s, lhs_violated_s,
	infected, cluster);
	for (int i = 0; i < rhs_violated_s.size(); i ++) {
		constraints_lhs_res.push_back(lhs_violated_s[i]);
		constraints_rhs_res.push_back(rhs_violated_s[i]);
	}

}

bool Separation::finds_constraints_COMBINATIONS(vector<double> infected) {
	//cout << "ALL COMBINATIONS SEPARATION SELECTED" << endl;
	//cout << "Begin separation" << endl;
	constraints_lhs_res.clear();
	constraints_rhs_res.clear();

	if (clusters.empty()) {
		cluster_find_all(infected);
	} else {
		for (int i = 0; i < clusters.size(); i ++) {
			cluster_find_all(infected, clusters[i]);
		}
	}
	//cout << "End separation" << endl;
	return ! constraints_rhs_res.empty();
}

bool Separation::finds_constraints_HEURISTIC(vector<double> infected) {
	//cout << "HEURISTIC SEPARATION SELECTED" << endl;
	//cout << "Begin separation" << endl;
	constraints_lhs_res.clear();
	constraints_rhs_res.clear();

	if (clusters.empty()) {
		initial_heuristic(infected);
	} else {
		for (int i = 0; i < clusters.size(); i ++) {
			initial_heuristic(infected, clusters[i]);
		}
	}
	//cout << "End separation" << endl;
	return ! constraints_rhs_res.empty();
}


bool Separation::finds_constraints(vector<double> infected) {
	return finds_constraints_HEURISTIC(infected);
}
