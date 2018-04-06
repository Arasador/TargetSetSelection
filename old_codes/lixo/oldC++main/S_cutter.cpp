#include "S_cutter.h"
#include <assert.h>

#typedef vector<int>::iterator Vector_it;

S_cutter::S_cutter(Vector2_int _adjacency_list, vector<int> _f) {
	adjacency_list = _adjacency_list;
	f = _f;
	N = f.size();
}

~S_cutter::S_cutter() {}

vector<int> S_cutter::infect_graph(vector<bool> infected) {
	vector<int> new_f(N, 0);
	vector<int> neighboors_infected(N, 0);
    vector<int> queue;
    for (int v = 0; v < N; v ++)
    	if (infected[v])
    		queue.push_back();
    while (! queue.empty()){
    	int v = queue[0];
    	queue.pop_front();
    	new_f[v] = -2;
    	for (Vector_it u = (adjacency_list[v]).begin(); u != (adjacency_list[v]).end(); ++u) {
    		if (infected[* u])
    			continue;
    		neighboors_infected[* u] ++;
    		if (-- new_f[* u] <= 0) {
    			queue.push_back(* u);
    			infected[* u] = True;
    		}
    	}
    }
    return new_f;
}

// finds | N(v) intersect (V/S)|
int S_cutter::neighbors_outside_s(vector<int> S, int v) {
	int solution = 0;
	vector<bool> S_b(N, false);
	for (int i = 0; i < S.size(); i ++)
		S_b[S[i]] = true;
	for (Vector_it u = (adjacency_list[v]).begin(); u != (adjacency_list[v]).end(); ++u) {
		if (! S_b[* u])
			solution ++;
	}
	return solution;
}


// for a given S, tries to find violated inequalities
int S_cutter::s_induced_inequality(vector<int> S) {
	assert(! S.empty());
	int min_S_infected, min_vertex = f[S[0]] + 1, -1;
    // finds min(f(v) - |N_G(v)) inter (V/S)|
	for (Vector_it v = S.begin(); v != S.end(); ++ v) {
		int new_min_value = f[v] - neighbors_outside_s(S, v);
		if (new_min_value < min_S_infected) {
			min_S_infected = new_min_value;
			min_vertex = v;
		}
	}
	//a vertex should be choosen as min
    assert(min_vertex != -1);
    return min_S_infected;
}


//given a set of infected vertices and an f, creates constraints based on
// the S set of remaining vertices
bool S_cutter::finds_S_constraints(vector<bool> infected, vector<int> new_f) {
	vector<bool> visited(N);
	for (int v = 0; v < N; v ++)
		visited[v] = new_f[v] <= 0;
	vector<int[2]>* components_limits = new vector<int[2]>();
	for (int w = 0; w < N; w ++) {
		if (visited[w])
			continue;
		vector<int> queue(1, w);
		visited[w] = true;
		vector<int> S;
		while (! queue.empty()) {
			int v = queue[0];
			queue.pop_front();
			if (new_f[v] > 0) {
                S.push_back(v)
                for (Vector_it u = adjacency_list[v].begin(); u != adjacency_list[v]; ++ u) {
                	if (! visited[* u]) {
                		queue.push_back(* u);
                		visited[* u] = true;
                	}
                }
			}
		}
		if (! S.empty()) {
            int S_lower_bound = S_induced_inequality(S);
            if (S_lower_bound > 0) {
                //assignment_constraint = cplex.SparsePair (ind = S,
                //    val = [1] * len(S))
                components_limits.push_back({S, S_lower_bound})
            }
		}
	}
	constraints = components_limits;
    return (! components_limits.empty())
}
