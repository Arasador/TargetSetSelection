/**
  S_cutter.h
  Given an S subset of of V set of vertices in a graph, we need to find 
  a new valid restriction, cutting the solution space, 
  If no such cut can be found, the optimal solution was found
**/
#ifndef SCUTTER
#define SCUTTER

#include "includes.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
#include <deque>
#include <assert.h>

using namespace std;

/*enum model {S_MODEL, S_SMALLER, WS_SMALLER, DOMINATED, WDOMINATED,
  S_SMALLER_H1, S_SMALLER_H2, S_SMALLER_NEW};
enum removals { REMOVED = -1, INFECTED = -2, DEGREE_ONE_REMOVAL = -3, FREE = 0};
//*/

//#define PRINT_LOG

class S_cutter {
  public:
  // variables
  vector<int> variables_used; // number of times total a variable was used
  vector<float> v_weights; // portion of total variable was used (sums 1)
  int variables_use_counter; // total a variables counter
  vector<vector<int> > adjacency_list;
  vector<vector<bool> > adjacency_matrix;
  vector<int> f, w;
  vector<bool> current_ub;
  int N, ub_val;

  // ----------------methods---------------

  int neighbors_outside_s(const vector<int> &S, int v);
  int S_induced_inequality(const vector<int> &S);
  bool finds_S_constraints(const vector<int> &new_f,
    vector<bool> visited = vector<bool>(), vector<double> weights = vector<double>());
  bool finds_s_model_constraints(vector<bool> &infected_vertices);
  // ----heuristic
  int heuristic_fmin_select_vertex(const vector<int> &new_f);
  int heuristic_degreefmin_select_vertex(const vector<int> &new_f);
  // -----weight
  void shuffle_positions_in_range(int i, int j, vector<int> &weights,
    vector<int> &positions);
  void shuffle_equal_positions(vector<int> &weights, vector<int> &positions);
  vector<int> weighted_positions(vector<int>&);
  vector<int> normal_shuffle();
  int select_next_vertex(const vector<int> &vertices,
    const vector<bool> &infected, int &position, model model_chosen = S_MODEL,
    vector<int> new_f = vector<int>());
  void reweight_vector_vertices_selected(const vector<vector<int> > &lhs);
  vector<int> weighted_option_selected(model model_chosen) ;
  // ----Smaller S
  int select_random_vertex(vector<bool> infected);
  void infect_one_vertex(int vi, vector<bool> &infected, vector<int> &new_f);
  bool can_select_random_component(vector<bool> infected, int selected_v);
    void shuffle_positions_in_range(int i, int j, vector<int> &weights);

  // Domination
  bool constraint_contained_in_other(const vector<bool> &constraint,
    const vector<bool> &bigger_constraint);
  bool constraint_is_dominated(const vector<int> &c_lhs_int, int c_rhs,
    const vector<int> &bc_lhs_int, int bc_rhs);
  void add_not_dominated_new_constraints(vector<vector<int> > &list_lhs,
    vector<int> &list_rhs);
  bool finds_s_with_domination_constraints(vector<bool> &infected, model model_chosen);
  // new S_smaller
  bool find_S_smaller_new_constraints(vector<bool> &infected, model model_chosen,
    vector<int> order = vector<int>(),
   vector<double> weights = vector<double>() );

  bool S_constraints_recursively(vector<bool> infected, vector<bool>&
    vertices_selected, vector<int> new_f, const vector<int>& selection_order,
    int position, vector<vector<int> >& found_constr_lhs,
    vector<int>& found_constr_rhs, int max_prev_rhs, bool& first, vector<bool>& upper_bound,
    vector<double>& weights);

  bool finds_all_components_S_small(vector<bool>& infected, model model_chosen, 
    vector<int> selection_order = vector<int>(), vector<double> weights = vector<double>());
  bool viable_new_right_side_constraint(vector<int>& max_rhs_constraints, int prev_rhs);
  // --- variables
  vector<vector<int> > constraints_lhs_res;
  vector<int> constraints_rhs_res;
  // --- constructor
  S_cutter(vector<vector<int> >, vector<int>, vector<int>);
  ~S_cutter();
  // methods
  bool finds_s_smaller_constraints(vector<bool> infected, model model_chosen);
  vector<int> infect_graph(vector<bool> &infected);
  bool finds_constraints(vector<bool> infected, model model_chosen);

  int initial_heuristic();
  int greedy_heuristic();


  // new function finds upper bounds
  void finds_upper_bound(int& upper_bound_weight, vector<bool>& ub_vertices,
    vector<bool>& infected, vector<int>& new_f, vector<int>& selection_order, 
    int& position, model model_chosen); 
};
#endif