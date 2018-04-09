#include "includes.h"

class S_cutter {
  public:
  // variables
  vector<int> variables_used; // number of times total a variable was used
  vector<float> v_weights; // portion of total variable was used (sums 1)
  int variables_use_counter; // total a variables counter
  vector<vector<int> > adjacency_list;
  vector<vector<bool> > adjacency_matrix;
  vector<int> f;
  int N;

  // ----------------methods---------------

  int neighbors_outside_s(const vector<int> &S, int v);
  int S_induced_inequality(const vector<int> &S);
  bool finds_S_constraints(const vector<int> &new_f);
  bool finds_s_model_constraints(vector<bool> &infected_vertices);
  // ----heuristic
  int heuristic_fmin_select_vertex(const vector<int> &new_f);
  int heuristic_degreefmin_select_vertex(const vector<int> &new_f);
  // -----weight
  void shuffle_positions_in_range(int i, int j, vector<int> &weights,
    vector<int> &positions);
  void shuffle_equal_positions(vector<int> &weights, vector<int> &positions);
  vector<int> weighted_positions();
  vector<int> normal_shuffle();
  int select_next_vertex(const vector<int> &vertices,
    const vector<bool> &infected, int &position, int option, vector<int> new_f);
  void reweight_vector_vertices_selected(const vector<vector<int> > &lhs);
  vector<int> weighted_option_selected(int option) ;
  // ----Smaller S
  int select_random_vertex(vector<bool> infected);
  void infect_one_vertex(int vi, vector<bool> &infected, vector<int> &new_f);
  bool can_select_random_component(vector<bool> infected);
    void shuffle_positions_in_range(int i, int j, vector<int> &weights);

  // Domination
  bool constraint_contained_in_other(const vector<bool> &constraint,
    const vector<bool> &bigger_constraint);
  bool constraint_is_dominated(const vector<int> &c_lhs_int, int c_rhs,
    const vector<int> &bc_lhs_int, int bc_rhs);
  void add_not_dominated_new_constraints(vector<vector<int> > &list_lhs,
    vector<int> &list_rhs);
  bool finds_s_with_domination_constraints(vector<bool> infected,  int option);
  void assert_valid_solution_found(vector<bool> infected);

  // --- variables
  vector<vector<int> > constraints_lhs_res;
  vector<int> constraints_rhs_res;
  // --- constructor
  S_cutter(Vector2_int, vector<int>);
  ~S_cutter();
  // methods
  bool finds_s_smaller_constraints(vector<bool> infected, int option);
  vector<int> infect_graph(vector<bool> &infected);
  bool finds_constraints(vector<bool> infected, int option);
};
//##############################################################################
//----------------- Auxiliar Functions

template<class T>
void print_vector(const vector<T> &v, string message) {
  cout << message;
  for (int i = 0; i < v.size(); i ++)
    cout << v[i] << "  ";
  cout << endl;
}

template<class T>
void print_matrix(const vector<vector<T> > &m, string message) {
  cout << message;
  for (int i = 0; i < m. size(); i ++)
    print_vector(m[i], "");
}

// detects if the option selected is a weighted option
bool is_weighted_option(int option) {
  return option == WS_SMALLER || option == WDOMINATED;
}
//##############################################################################
// ----------------Constructor 

S_cutter::S_cutter(Vector2_int _adjacency_list, vector<int> _f) {
  adjacency_list = _adjacency_list;
  f = _f;
  N = f.size();
  variables_used = vector<int>(N);
  v_weights = vector<float>(N);
  variables_use_counter = 0;
}

S_cutter::~S_cutter() {}

//##############################################################################
// ----------------Functions related to S

// given a vector with all the infected vertices, spreads the infection and
// returns a new f vector, where vertices infected have value INFECTED
vector<int> S_cutter::infect_graph(vector<bool> &infected) {
  //copies f to new_f, to start creating the new f
  vector<int> new_f(f.begin(), f.end());
  deque<int> queue;
  //put already infected vertices in tue queue, so they can spread the infection 
  for (int v = 0; v < N; v ++) {
    if (infected[v]) {
      queue.push_back(v);
      new_f[v] = INFECTED;
    }
  }
  while (! queue.empty()) {
    int v = queue[0];
    queue.pop_front();
  // for each infected vertex, see if neighbors are infected 
    for (vector<int>::iterator u = (adjacency_list[v]).begin();
      u != (adjacency_list[v]).end(); ++ u) {
      if (infected[* u])
        continue;
      // infects new vertex, puts it in the queue
      if ((-- new_f[* u]) <= 0) {
        infected[* u] = true;
        queue.push_back(* u);
        new_f[* u] = INFECTED;
      }
    }
  }
  return new_f;
}

// finds | N(v) intersect (V/S)|
int S_cutter::neighbors_outside_s(const vector<int> &S, int v) {
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

// for a given S, tries to find violated inequalities
int S_cutter::S_induced_inequality(const vector<int> &S) {
  assert(! S.empty());
  int S_lower_bound = f[S[0]] - neighbors_outside_s(S, S[0]);
  for (int v = 1; v < S.size(); v ++) {
    int current_value = f[S[v]] - neighbors_outside_s(S, S[v]);
    if (current_value < S_lower_bound)
      S_lower_bound = current_value;
  }
  return S_lower_bound;
}

//given a set of infected vertices and an f, creates constraints based on
// the S set of remaining vertices
bool S_cutter::finds_S_constraints(const vector<int> &new_f) {
  // marks infected vertices as visited
  vector<bool> visited(N);
  for (int v = 0; v < N; v ++)
    visited[v] = (new_f[v] <= 0);
  
  vector<int> rhs;
  vector<vector<int> > lhs;
  for (int w = 0; w < N; w ++) {
    if (visited[w])
      continue;
    // if selects someone, a new S has to be tested for constraints
    vector<bool> selected_in_S(N, false);
    bool S_is_not_empty = false;
    deque<int> queue;
    queue.push_back(w);
    visited[w] = true;
    // for each w not visited,finds its connected component
    while(! queue.empty()) {
      int v = queue[0];
      queue.pop_front();
      // only adds vertices active
      if (new_f[v] <= 0)
        continue;
      selected_in_S[v] = true;
      S_is_not_empty = true;
      // breath first search spreading to neighbors not visited 
      //(which should be all)
      for (vector<int>::iterator u = (adjacency_list[v]).begin();
        u != (adjacency_list[v]).end(); ++ u) {
        if (visited[* u] || new_f[* u] <= 0)
          continue;
        queue.push_back(* u);
        visited[* u] = true;
      }
    }
    // given an S, finds new constraint 
    if (S_is_not_empty) {
      vector<int> S;
      for (int i = 0; i < N; i ++)
        if (selected_in_S[i])
          S.push_back(i);

      int S_lower_bound = S_induced_inequality(S);
      if (S_lower_bound > 0) {
        lhs.push_back(S);
        rhs.push_back(S_lower_bound);
      }
    }
  }

  if (rhs.empty()) {
    return false;
  }
  constraints_rhs_res = rhs;
  constraints_lhs_res = lhs;
  return true;
}
//##############################################################################
//----------- WEIGHTED METHODS

// given a range [i,j), shuffle its positions (used so equally ranked positions
//dont always end up in the same order)
void S_cutter::shuffle_positions_in_range(int i, int j, vector<int> &weights,
  vector<int> &positions){
  vector<int> new_positions(positions.begin() + i, positions.begin() + j);
  random_shuffle(new_positions.begin(), new_positions.end());
  for (int n = i; n < j; n ++)
    positions[n] = new_positions[n - i];
}

// finds each range of positions in the same rank, and shuffles them
void S_cutter::shuffle_equal_positions(vector<int> &weights, vector<int> 
	&positions) {
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
    assert(weights[i - 1] >= weights[i]);
}

// auxiliary function for the sorting order we want, when ordering pairs
bool sort_pairs (pair<int, int> i,pair<int, int> j) {
  return (i.second > j.second);
}

//finds a order of vertices to be selected depending on their "weight" (number
// of times they were already used in constraints)
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

// shuffles vertices randomly
vector<int> S_cutter::normal_shuffle(){
  vector<int> v(N);
  for (int i = 0; i < N; i ++)
    v[i] = i;
  random_shuffle(v.begin(), v.end());
  return v;
}

// given an order of vertices, gives next that is not infected
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
void S_cutter::reweight_vector_vertices_selected(const vector<vector<int> > 
	&lhs) {
  for (int i = 0; i < lhs.size(); i ++) {
    for (int j = 0; j < lhs[i].size(); j ++) {
      variables_used[lhs[i][j]] ++;
      variables_use_counter ++;
    }
  }
  // if counters gets to big, reduces it (but we will have bigger problems)

  if (variables_use_counter > 100000) {
    for (int i = 0; i < N; i ++)
      variables_used[i] = variables_used[i] / 1000;
  }
  if (variables_use_counter == 0)
    return;
  // updates v_weights
  for (int i = 0; i < N; i ++)
    v_weights[i] = float(variables_used[i]) / variables_use_counter;
}

// sees which weighted option is selected
vector<int> S_cutter::weighted_option_selected(int option) {
  if (option == WS_SMALLER || option == WDOMINATED) {
    //return weighted_shuffle()
    return weighted_positions();
  }
  return normal_shuffle();
}
//##############################################################################
//------ FIRST SIMPLE S MODEL

// gets the not infected vertices, uses that as s and creates constraints
bool S_cutter::finds_s_model_constraints(vector<bool> &infected_vertices){
  vector<int> new_f = infect_graph(infected_vertices);
  return finds_S_constraints(new_f);
}
//##############################################################################
//# ------ SMALLER S
//# NEW CUTTING ALGORITHM THAT INFECTS THE REST and gets smaller S

// returns next vertex to be infected, that is a neighbor of a vertex
// that has the least using heuristic
int S_cutter::heuristic_fmin_select_vertex(const vector<int> &new_f) {
  int min_f = INFECTED, min_arg = -1;

  for (int v = 0; v < N; v ++) {
    if (new_f[v] < 0)
      continue;
    if (new_f[v] < min_f || min_f == INFECTED){
      min_f = new_f[v];
      min_arg = v;
    }
  }
  //cout << min_arg << endl;
  assert(min_arg >= 0);
  return min_arg;
}

// returns next vertex to be infected, that is a neighbor of a vertex
// that has the least using heuristic
int S_cutter::heuristic_degreefmin_select_vertex(const vector<int> &new_f) {
  float min_f = INFECTED;
  int min_arg = -1;
  for (int v = 0; v < N; v ++) {
    if (new_f[v] == INFECTED)
      continue;
    float aux = (float) new_f[v] / (float) adjacency_list.size();
    if (aux < min_f || min_f == INFECTED) {
      min_f = new_f[v];
      min_arg = v;
    }
  }
  assert(min_arg >= 0);
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
  new_f[vi] = INFECTED;
  while (! queue.empty()) {
    int v = queue[0];
    queue.pop_front();
    for (Vector_it u = adjacency_list[v].begin(); u != adjacency_list[v].end();
      ++ u) {
      if (infected[* u])
        continue;
      if ((-- new_f[* u]) <= 0) {
        queue.push_back(* u);
        infected[* u] = true;
        new_f[* u] = INFECTED;
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
  vector<int> new_f = infect_graph(infected);
  if (! finds_S_constraints(new_f))
  	return false;
  // gets a vector with the section order of vertices based on the option
  // weighted or not
  vector<int> selection_order = weighted_option_selected(option);
  //print(selection_order)
  int position = 0;
  //keeps selecting one of the avaliable components
  while (can_select_random_component(infected)) {
    int v = select_next_vertex(selection_order, infected, position, option, 
    	new_f);
    assert(v != -1);
    infect_one_vertex(v, infected, new_f);
    if (! finds_S_constraints(new_f)){
      break;
    }
  }
  if (is_weighted_option(option)) {
    reweight_vector_vertices_selected(constraints_lhs_res);
  }
  for (int r = 0; r < constraints_rhs_res.size(); r ++)
    assert(constraints_rhs_res[r] == 1);
  return constraints_rhs_res.size() != 0;
}
//##############################################################################
//# --------- DOMINATION
// tests if given 2 constraint, one is constained by the others
bool S_cutter::constraint_contained_in_other(const vector<bool> &constraint,
  const vector<bool> &bigger_constraint) {
  assert(constraint.size() == N);
  for (int i = 0; i < N; i ++)
    if (constraint[i] and ! bigger_constraint[i])
      return false;
  return true;
}

// returns if a given constraint already is represented, by being dominated by
// another constraint
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
      left_difference += 1 - c_lhs[i]; 
    // difference in right side of constraints
  }
  int right_difference = bc_rhs - c_rhs;
  // if there are to many variables removed, compared with the min number
  // of activated vertices mandatory, then it is not dominated and must be
  // added this new constraint
  return left_difference <= right_difference;
}

// addes not dominated constraints
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
//*/

//this function returns if we found constraints for a smaller S and domination
bool S_cutter::finds_s_with_domination_constraints(vector<bool> infected,
  int option) {
  vector<int> new_f = infect_graph(infected);
  if (!finds_S_constraints(new_f))
  	return false;
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
    if (finds_S_constraints(new_f))
      add_not_dominated_new_constraints(list_lhs, list_rhs);
    else
      break;
    if (is_weighted_option(option)) {
      reweight_vector_vertices_selected(list_lhs);
      //print("v_vector ", self.variables_used )
    }
  }
  constraints_lhs_res = list_lhs;
  constraints_rhs_res = list_rhs;
  return ! constraints_lhs_res.empty();
}

// 
void S_cutter::assert_valid_solution_found(vector<bool> infected) {
  vector<int> new_f = infect_graph(infected);
  for (int i = 0;  i < N; i ++)
    assert(new_f[i] == INFECTED);
}


//##############################################################################
// selectes which model will run
bool S_cutter::finds_constraints(vector<bool> infected, int option) {

  switch (option) {
    case S_MODEL:
      //cout << "---------------S_model-----------------" << endl;
      return finds_s_model_constraints(infected);
    case S_SMALLER:
      //cout << "---------------S_smaller-----------------" << endl;
      return finds_s_smaller_constraints(infected, option);
    case WS_SMALLER:
      //cout << "---------------S_smallerWW-----------------" << endl;
      return finds_s_smaller_constraints(infected, option);
    case S_SMALLER_H1:
      //cout << "here H1" << endl;
      return finds_s_smaller_constraints(infected, option);
    case S_SMALLER_H2:
      //cout << "here H2" << endl;
    return finds_s_smaller_constraints(infected, option);
    case DOMINATED:
      //cout << "---------------Dominated-----------------" << endl;
      return finds_s_with_domination_constraints(infected, option);
    case WDOMINATED:
      //cout << "---------------DominatedWWWW-----------------" << endl;
      return finds_s_with_domination_constraints(infected, option);
      //*/
  }
}
