#ifndef SEPARATION
#define SEPARATION

#include <iostream>
#include <algorithm>
#include <vector>


using namespace std;

class Separation {
public:
	vector<vector<int>> adjacency_list, constraints_lhs_res;
	vector<vector<bool>> clusters;
	vector<int> f, w, num_neightbors, constraints_rhs_res;
	vector<double> costs;
	int N;
	Separation(vector<vector<int>> _adjacency_list, vector<int> _f,
		vector<int> _w);
	~Separation();
	int neighbors_outside_s(const vector<int> &S, int v);
	int neighbors_outside_s_bool(const vector<bool> &S, int v);
	int S_induced_inequality(const vector<int> &S);
	int S_induced_inequality_bool(const vector<bool> &S);
	int initial_vertex_selection(vector<double>& infected,
		 const vector<bool>& cluster = vector<bool>());
	void add_to_neighbors(vector<bool>& neighbors_S, int v0,
		const vector<bool>& cluster);
	int get_max_S_neighbor(const vector<bool> &neighbors, vector<bool>& S,
		int& rhs, const vector<double>& infected);
	void initial_heuristic(vector<double>& infected,
		const vector<bool>& cluster = vector<bool>());
	void check_if_violates(vector<int> S, vector<int>& rhs,
		vector<vector<int>>& lhs, vector<double>& infected);
	void generates_Ss_and_tests(int i, vector<int>& S, int N, vector<int>& rhs,
		vector<vector<int>>& lhs, vector<double>& infected,
		const vector<bool>& cluster = vector<bool>());
	void tests_all_Ss(vector<double>& infected);
	void cluster_find_all(vector<double>& infected,
		const vector<bool>& cluster = vector<bool>());
	bool finds_constraints(vector<double> infected);
	bool finds_constraints_COMBINATIONS(vector<double> infected);
	bool finds_constraints_HEURISTIC(vector<double> infected);
};
#endif