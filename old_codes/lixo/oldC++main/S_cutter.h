//#ifndef SCUTTER_H
//#define SCUTTER_H
#include<assert.h>
#include <iostream>
#include <vector>
#include<deque>
#include <map>
#include <stack>
#include<string>
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

typedef vector<vector<int> > Vector2_int;
typedef vector<vector<bool> > Vector2_bool;
typedef vector<int>::iterator Vector_it;

class S_cutter {
protected:
	vector<int> variables_used;
	vector<float> v_weights;
	int counter;

public:
	vector<vector<int> > adjacency_list, constraints_lhs_res;
	vector<vector<bool> > adjacency_matrix;
	vector<int> f, constraints_rhs_res;
	int N;
	S_cutter(Vector2_int, vector<int>);
	~S_cutter();
	vector<int> infect_graph(vector<bool> infected);
	int neighbors_outside_s(vector<int> S, int v);
	int S_induced_inequality(vector<int> S);
	bool finds_S_constraints(const vector<bool> &infected, vector<int> &new_f);
	bool finds_s_model_constraints(vector<bool> &infected_vertices);
	bool finds_constraints(vector<bool> infected, int option);
	// heuristic
	int heuristic_fmin_select_vertex(const vector<int> &new_f);
	int heuristic_degreefmin_select_vertex(const vector<int> &new_f);
	// weight
	void shuffle_positions_in_range(int i, int j, vector<int> &weights,
		vector<int> &positions);
	void shuffle_equal_positions(vector<int> &weights, vector<int> &positions);
	vector<int> weighted_positions();
	vector<int> normal_shuffle();
	int select_next_vertex(const vector<int> &vertices,
		const vector<bool> &infected, int &position, int option, vector<int> new_f) ;
	void reweight_vector_vertices_selected(const vector<vector<int> > &lhs);
	vector<int> weighted_option_selected(int option) ;

	// Smaller S
	int select_random_vertex(vector<bool> infected);
	void infect_one_vertex(int vi, vector<bool> &infected, vector<int> &new_f);
	bool can_select_random_component(vector<bool> infected);
	bool finds_s_smaller_constraints(vector<bool> infected, int option);
	void shuffle_positions_in_range(int i, int j, vector<int> &weights);

	// Domination
	bool constraint_contained_in_other(const vector<bool> &constraint,
		const vector<bool> &bigger_constraint);
	bool constraint_is_dominated(const vector<int> &c_lhs_int, int c_rhs,
		const vector<int> &bc_lhs_int, int bc_rhs);
	void add_not_dominated_new_constraints(vector<vector<int> > &list_lhs,
		vector<int> &list_rhs);
	bool finds_s_with_domination_constraints(vector<bool> infected,	int option);
};

template<class T>
void print_vector(vector<T> &v, string message) {
  cout << message;
  for (int i = 0; i < v.size(); i ++)
    cout << v[i] << "  ";
  cout << endl;
}

template<class T>
void print_matrix(vector<vector<T> > &m, string message) {
  cout << message;
  for (int i = 0; i < m. size(); i ++)
    print_vector(m[i], "");
}

bool is_weighted_option(int option) {
	return option == WS_SMALLER || option == WDOMINATED;
}

S_cutter::S_cutter(Vector2_int _adjacency_list, vector<int> _f) {
	adjacency_list = _adjacency_list;
	f = _f;
	N = f.size();
	variables_used = vector<int>(N);
	v_weights = vector<float>(N);
	counter = 0;
}

S_cutter::~S_cutter() {}

vector<int> S_cutter::infect_graph(vector<bool> infected) {
	vector<int> new_f(f.begin(), f.end());
	vector<int> neighboors_infected(N, 0);
    deque<int> queue;
    for (int v = 0; v < N; v ++)
    	if (infected[v])
    		queue.push_back(v);

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
    			infected[* u] = true;
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
int S_cutter::S_induced_inequality(vector<int> S) {
	assert(! S.empty());
	int min_S_infected = f[S[0]] + 1;
	int min_vertex = -1;
    // finds min(f(v) - |N_G(v)) inter (V/S)|
	for (Vector_it v = S.begin(); v != S.end(); ++ v) {
		int new_min_value = f[* v] - neighbors_outside_s(S, * v);
		if (new_min_value < min_S_infected) {
			min_S_infected = new_min_value;
			min_vertex = * v;
		}
	}
	//a vertex should be choosen as min
    assert(min_vertex != -1);
    return min_S_infected;
}


//given a set of infected vertices and an f, creates constraints based on
// the S set of remaining vertices
bool S_cutter::finds_S_constraints(const vector<bool> &infected,
	vector<int> &new_f) {
	vector<bool> visited(N);
	constraints_lhs_res = vector<vector<int> >();
	constraints_rhs_res = vector<int>();
	for (int v = 0; v < N; v ++)
		visited[v] = new_f[v] <= 0;
	for (int w = 0; w < N; w ++) {
		if (visited[w])
			continue;
		deque<int> queue(1, w);
		visited[w] = true;
		vector<int> S;
		while (! queue.empty()) {
			int v = queue[0];
			queue.pop_front();
			if (new_f[v] > 0) {
        S.push_back(v);
        for (Vector_it u = adjacency_list[v].begin(); u != (adjacency_list[v]).end(); ++ u) {
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
                constraints_lhs_res.push_back(S);
								constraints_rhs_res.push_back(S_lower_bound);
            }
		}
	}
    return (! constraints_rhs_res.empty());
}

//----------- WEIGHTED METHODS



void S_cutter::shuffle_positions_in_range(int i, int j, vector<int> &weights,
	vector<int> &positions){
	vector<int> new_positions(positions.begin() + i, positions.begin() + j);
	random_shuffle(new_positions.begin(), new_positions.end());
	for (int n = i; n < j; n ++)
		positions[n] = new_positions[n - i];
}

void S_cutter::shuffle_equal_positions(vector<int> &weights, vector<int> &positions) {
	int current_weight = weights[0];
	int i = 0;
	for (int j = 1; j < N; j ++) {
		if (current_weight != weights[j]) {
			shuffle_positions_in_range(i, j, weights, positions);
			i = j;
			current_weight = weights[j];
		}
	}
	shuffle_positions_in_range(i, N, weights, positions);
	//asserts everything is in order
	for (int i = 1; i < N; i ++)
		assert(weights_vector[i - 1][1] >= weights_vector[i][1]);
}

bool sort_pairs (pair<int, int> i,pair<int, int> j) {
	return (i.second > j.second);
}

vector<int> S_cutter::weighted_positions() {
	vector<pair<int, int> > pos_weights(N);
	for (int i = 0; i < N; i ++) {
		pos_weights[i] = make_pair(i, variables_used[i]);
	}
	sort(pos_weights.begin(), pos_weights.end(), sort_pairs);

	vector<int> positions(N), weights(N);;
	for (int i = 0; i < N; i ++) {
		positions[i] = pos_weights[i].first;
		weights[i] = pos_weights[i].second;
	}
	shuffle_equal_positions(weights, positions);
	return positions;
}


vector<int> S_cutter::normal_shuffle(){
	vector<int> r(N);
	for (int i = 0; i < N; i ++)
		r[i] = i;
	random_shuffle(r.begin(), r.end());
	return r;
}


int S_cutter::select_next_vertex(const vector<int> &vertices,
	const vector<bool> &infected, int &position, int option, vector<int> new_f) {
	if (option == S_SMALLER_H1) {
		return heuristic_fmin_select_vertex(new_f);
	}
	else if (option == S_SMALLER_H2) {
		return heuristic_degreefmin_select_vertex(new_f);
	}

	for (int i = position; i < vertices.size(); i ++) {
		if (!infected[vertices[i]]) {
			position = i;
			return vertices[i];
		}
	}
	return -1;
}

// When new constraints are found, reweight the vector that counts variables
// used in the constraints
void S_cutter::reweight_vector_vertices_selected(const vector<vector<int> > &lhs) {
	for (int i = 0; i < lhs.size(); i ++) {
		for (int j = 0; j < lhs[i].size(); j ++) {
			variables_used[lhs[i][j]] ++;
			counter ++;
		}
	}
	// if counters gets to big, reduces it (but we will have bigger problems)

	if (counter > 100000) {
		for (int i = 0; i < N; i ++)
			variables_used[i] = variables_used[i] / 1000;
	}
	if (counter == 0)
		return;
	// updates v_weights
	for (int i = 0; i < N; i ++)
		v_weights[i] = float(variables_used[i]) / counter;
}


vector<int> S_cutter::weighted_option_selected(int option) {
	if (option == WS_SMALLER || option == WDOMINATED) {
		//return weighted_shuffle()
		return weighted_positions();
	}
	return normal_shuffle();
}


//------ FIRST SIMPLE S MODEL

// gets the not infected vertices, uses that as s and creates constraints
bool S_cutter::finds_s_model_constraints(vector<bool> &infected_vertices){
	vector<int> new_f = infect_graph(infected_vertices);
	//print_vector(new_f, "New f found");
	return finds_S_constraints(infected_vertices, new_f);
}

//##############################################################################
//# ------ SMALLER S
//# NEW CUTTING ALGORITHM THAT INFECTS THE REST and gets smaller S

// returns next vertex to be infected, that is a neighbor of a vertex
// that has the least using heuristic


int S_cutter::heuristic_fmin_select_vertex(const vector<int> &new_f) {
  int min_f = numeric_limits<int>::max(), min_arg = -1;
  for (int v = 0; v < N; v ++) {
  	if (new_f[v] == -2)
  	  continue;
  	if (new_f[v] < min_f){
  		min_f = new_f[v];
  		min_arg = v;
  	}
  }
  return min_arg;
}

int S_cutter::heuristic_degreefmin_select_vertex(const vector<int> &new_f) {
  float min_f = numeric_limits<float>::max();
  int min_arg = -1;
  for (int v = 0; v < N; v ++) {
  	if (new_f[v] == -2)
  	  continue;
  	float aux = (float) new_f[v] / (float) adjacency_list.size();
  	if (aux < min_f) {
  		min_f = new_f[v];
  		min_arg = v;
  	}
  }
  return min_arg;
}

// given a range and a set of forbidden vertices, selects the one avaliable
// if no vertex is selected, returns -1
int S_cutter::select_random_vertex(vector<bool> infected) {
	int num_not_infected = 0;
	for (int i = 0; i < N; i ++)
		if (! infected[i])
			num_not_infected ++;
	int r = rand() % num_not_infected;
	int j = 0;
	for (int i = 0; i < N; i ++) {
		if (infected[i])
			continue;
		if (j ++ == r)
			return i;
	}
	return -1;
}
// given one new infected vertex in the graph, carries on the infection
// this vertex cannot be already infected
void S_cutter::infect_one_vertex(int vi, vector<bool> &infected,
	vector<int> &new_f) {
  deque<int> queue;
	queue.push_back(vi);
  assert(! infected[vi]);
  infected[vi] = true;
  new_f[vi] = -2;
  while (! queue.empty()) {
		int v = queue[0];
		queue.pop_front();
		for (Vector_it u = adjacency_list[v].begin(); u != adjacency_list[v].end();
			++ u) {
			if (infected[* u])
				continue;
			if (-- new_f[* u] <= 0) {
				queue.push_back(* u);
				infected[* u] = true;
				new_f[* u] = -2;
			}
		}
	}
}

// modifies infected vertices, so everything outside one component counts as
// infected. Return if we still have components
bool S_cutter::can_select_random_component(vector<bool> infected) {
	int cv = -1;
	for (int i = 0; i < N; i ++)
		if (! infected[i]) {
			cv = i;
			break;
		}
	if (cv == -1)
		return false;
	deque<int> queue;
	queue.push_back(cv);
	vector<bool> component(N, false);
	while (! queue.empty()) {
		int v = queue[0];
		queue.pop_front();
		component[v] = true;
		for (Vector_it u = adjacency_list[v].begin(); u != adjacency_list[v].end();
		++ u) {
			if (! infected[* u] && !component[* u]) {
				component[* u] = true;
				queue.push_back(* u);
			}
		}
	}
	// modifies infected vertices, so if v is outside component, treat as
	// infect
	for (int i = 0; i < N; i ++)
		infected[i] = infected[i] || ! component[i];
	return true;
}



//this function returns if we found constraints for a smaller S
bool S_cutter::finds_s_smaller_constraints(vector<bool> infected, int option) {
	//cout << "inside scutter" << endl;
  vector<int> new_f = infect_graph(infected);
  print_vector(new_f, "new_f ");
	finds_S_constraints(infected, new_f);
  // gets a vector with the section order of vertices based on the option
  // weighted or not
  vector<int> selection_order = weighted_option_selected(option);
  //print(selection_order)
  int position = 0;
	vector<vector<int> > lhs;
	vector<int> rhs;
  //keeps selecting one of the avaliable components
  while (can_select_random_component(infected)) {
		cout << "entered" << endl;
      //v = self.select_random_vertex(already_infected)
      int v = select_next_vertex(selection_order, infected, position, option,
				new_f);
			cout << "v " << v << endl;
      assert(v != -1);
      infect_one_vertex(v, infected, new_f);
			print_vector(new_f, "new f infected after v ");
      if (finds_S_constraints(infected, new_f)){
				lhs = constraints_lhs_res;
				rhs = constraints_rhs_res;
			}
			cout << "found constraints" << endl;
		}
  if (is_weighted_option(option)) {
		cout << "reweight " << endl;
      reweight_vector_vertices_selected(lhs);
      //print("v_vector ", self.variables_used )
	}
	cout << "left " << endl;
  //print("Final bounds ", bounds)
	constraints_lhs_res = lhs;
  constraints_rhs_res = rhs;
	if (! lhs.empty())
		print_vector(lhs[0], "test");
	cout << "out" << endl;
  return lhs.size() != 0;
}



//##############################################################################
//# --------- DOMINATION
bool S_cutter::constraint_contained_in_other(const vector<bool> &constraint,
	const vector<bool> &bigger_constraint) {
	assert(constraint.size() == N);
	for (int i = 0; i < N; i ++)
		if (constraint[i] and ! bigger_constraint[i])
			return false;
	return true;
}

bool S_cutter::constraint_is_dominated(const vector<int> &c_lhs_int, int c_rhs,
	const vector<int> &bc_lhs_int, int bc_rhs) {
	// constraint lhs bool, bigger constraint lhs
	vector<bool> c_lhs(N, false), bc_lhs(N, false);
	// convert constraint in int to bool
	for (int i = 0; i < c_lhs_int.size(); i ++)
		c_lhs[c_lhs_int[i]] = true;
	for (int i = 0; i < bc_lhs_int.size(); i ++)
		bc_lhs[bc_lhs_int[i]] = true;
	// if constraint contained id bigger one, return is dominated
	if (! constraint_contained_in_other(c_lhs, bc_lhs))
		return false;

	int left_difference = 0;
	for (int i = 0; i < N; i ++) {
		if (bc_lhs[i])
			// number of different variables between each constraint
			left_difference += 1 - c_lhs[i]; // if have variable in bc, if c not has, add
		// difference in right side of constraints
	}
	int right_difference = bc_rhs - c_rhs;
	// if there are to many variables removed, compared with the min number
	// of activated vertices mandatory, then it is not dominated and must be
	// added this new constraint
	return left_difference <= right_difference;
}

void S_cutter::add_not_dominated_new_constraints(vector<vector<int> > &new_lhs,
	vector<int> &new_rhs) {
	int num_constraints = new_lhs.size();
	for (int i = 0; i < constraints_lhs_res.size(); i ++) {
		bool constraint_not_dominated = true;
		for (int j = 0; j < num_constraints; j ++) {
			if (constraint_is_dominated(constraints_lhs_res[i], constraints_rhs_res[i], new_lhs[i], new_rhs[i])) {
				constraint_not_dominated = false;
				break;
			}
		}
		if (constraint_not_dominated) {
			new_lhs.push_back(constraints_lhs_res[i]);
			new_rhs.push_back(constraints_rhs_res[i]);
		}
	}
}

/*void S_cutter::add_not_dominated_new_constraints_smaller(vector<vector<int> >
	&list_lhs, vector<int> &list_rhs) {
	int num_constraints = list_lhs.size();
	for (int j = 0; j < num_constraints; j ++) {
		bool constraint_not_dominated = true;
		for (int i = 0; i < lhs.size(); i ++) {
			if (constraint_is_dominated(lhs[i], rhs[i], list_lhs[i], list_rhs[i])) {
				constraint_not_dominated = false;
				break;
			}
		}
		if (constraint_not_dominated) {
			list_lhs.push_back(lhs[i]);
			list_rhs.push_back(rhs[i]);
		}
	}
}

        new_bounds = [_ for _ in self.constraints]
        #print(new_bounds)
        for bigger_constraint in bounds:
            constraint_not_dominated = True
            for constraint in self.constraints:
                if self.constraint_is_dominated(constraint, bigger_constraint):
                    constraint_not_dominated = False
                    break
            if constraint_not_dominated:
                new_bounds.append(bigger_constraint)
        bounds += new_bounds
//*/

//this function returns if we found constraints for a smaller S and domination
bool S_cutter::finds_s_with_domination_constraints(vector<bool> infected,
	int option) {
	vector<int> new_f = infect_graph(infected);
	finds_S_constraints(infected, new_f);
	vector<vector<int> > list_lhs = constraints_lhs_res;
	vector<int> list_rhs = constraints_rhs_res;
	// gets a vector with the section order of vertices based on the option
	// weighted or not
	vector<int> selection_order = weighted_option_selected(option);
	int position = 0;
	// keeps selecting one of the avaliable components
	while (can_select_random_component(infected)) {
		int v = select_next_vertex(selection_order, infected, position, option, new_f);
		assert(v != -1);
		infect_one_vertex(v, infected, new_f);
		if (finds_S_constraints(infected, new_f))
			add_not_dominated_new_constraints(list_lhs, list_rhs);

		if (is_weighted_option(option)) {
			reweight_vector_vertices_selected(list_lhs);
			//print("v_vector ", self.variables_used )
		}
	}
	constraints_lhs_res = list_lhs;
	constraints_rhs_res = list_rhs;
	return ! constraints_lhs_res.empty();
}





bool S_cutter::finds_constraints(vector<bool> infected, int option) {
	//cout << "try find constraint" << endl;
	//print_vector(infected, "Infected: ");

	switch (option) {
		case S_MODEL:
			return finds_s_model_constraints(infected);
		case S_SMALLER:
			return finds_s_smaller_constraints(infected, option);
    case WS_SMALLER:
			return finds_s_smaller_constraints(infected, option);
		case S_SMALLER_H1:
			cout << "here H1" << endl;
			return finds_s_smaller_constraints(infected, option);
		case S_SMALLER_H2:
		cout << "here H2" << endl;
		return finds_s_smaller_constraints(infected, option);
    case DOMINATED:
			return finds_s_with_domination_constraints(infected, option);
    case WDOMINATED:
			return finds_s_with_domination_constraints(infected, option);
			//*/
  }
}
