from __future__ import print_function

import sys
import time
import re
import math
from s_cutter import S_cutter

import cplex
from cplex.exceptions import CplexSolverError,CplexError
from cplex.callbacks import LazyConstraintCallback
from cplex.callbacks import UserCutCallback, LazyConstraintCallback



def assert_valid_solution(s_cutter, infected_vertices):
    new_f = s_cutter.infect_graph(infected_vertices)
    for i in new_f:
        assert(i < 0)
    print("asserted")

#second one, where you find smaller S, and uses this constraint
class LazySCallback(LazyConstraintCallback):
    def __call__(self):
        if self.has_incumbent():
            self.gap = min(self.gap, self.get_MIP_relative_gap())
        f = self.f
        s_cutter = self.s_cutter
        infected_vertices = [bool(self.get_values(v)) for v in range(len(f))]
        found_constraints = False

        if s_cutter.finds_constraints(infected_vertices, self.s_option):
            #self.tolerances.uppercutoff = min(s_cu)
            for constraint in s_cutter.constraints:
                self.constraints_counter += 1
                assignment_constraint = cplex.SparsePair (ind = constraint[0],
                     val = [1] * len(constraint[0]))
                S_lb = constraint[1]
                self.add(constraint=assignment_constraint,
                                             sense="G",
                                             rhs=S_lb)

        else:
            #print("Ended search tree")
            #print(infected_vertices)
            new_f = s_cutter.infect_graph(infected_vertices)
            for i in new_f:
                assert(i < 0)

        self.end_time = self.get_time()

#initial S with S = V
def add_initial_constraint(model, N, f, w):
    #initial infected state of each vertex
    x = range(N)
    #set objective function
    model.objective.set_sense(model.objective.sense.minimize)
    model.variables.add(obj = [1] * N,#w,
                        lb = [0] * N,
                        ub = [1] * N,
                        types = ["B"] * N)
    #and initial inequality, where sum(x) >= min(f)
    assignment_constraint = cplex.SparsePair (ind = x, val = [1] * N)
    model.linear_constraints.add(lin_expr=[assignment_constraint],
                                 senses=["G"],
                                 rhs=[min(f)])

def PCI_model_2(adjacencylist, f, w, s_option, timelimit):
    N = len(f)
    model = cplex.Cplex()
    add_initial_constraint(model, N, f, w)
    #create object for the S_cutter
    s_cutter = S_cutter(adjacencylist, f)

    model.parameters.preprocessing.presolve.set(
        model.parameters.preprocessing.presolve.values.off)

    #model.parameters.threads.set(1)

    model.parameters.mip.strategy.search.set(
        model.parameters.mip.strategy.search.values.traditional)

    #defines lazy callback and parameters it has
    lazy_s_cb = model.register_callback(LazySCallback)
    lazy_s_cb.adjacencylist = adjacencylist
    lazy_s_cb.f = f
    lazy_s_cb.s_cutter = s_cutter
    lazy_s_cb.s_option = s_option
    lazy_s_cb.starttime = model.get_time()
    # starts with one constraint
    lazy_s_cb.constraints_counter = 1
    lazy_s_cb.gap = 1

    model.parameters.timelimit.set(timelimit)
    model.parameters.preprocessing.presolve.set(1)
    #model.parameters.parallel.set(-1)


    #model.parameters.parallel
    #model.set_results_stream(None)
    start_time = model.get_time()
    try:
        model.solve()
    except CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        solution = model.solution
        sol_value = solution.get_objective_value()
        #gap = solution.get_mip_relative_gap()
        print("GAAP: ", lazy_s_cb.gap)
        #time = float('%.4f'%(lazy_s_cb.end_time - start_time))
        time = model.get_time() - start_time
        print(sol_value)
        infected_vertices = [bool(solution.get_values(v)) for v in range(N)]
        assert_valid_solution(s_cutter, infected_vertices)
        print("\nSolution status = ", solution.get_status(), ":", end=' ')
        print(solution.status[solution.get_status()])

        return [sol_value, time, lazy_s_cb.constraints_counter,
            solution.get_status(), lazy_s_cb.gap]

def s_model_initial(adjacencylist, f, w, timelimit):
    return PCI_model_2(adjacencylist, f, w, "s_model", timelimit)

def smallest_s_model(adjacencylist, f, w, timelimit):
    return PCI_model_2(adjacencylist, f, w, "smaller_s", timelimit)

# w is for weighted
def smallest_s_weighted(adjacencylist, f, w, timelimit):
    return PCI_model_2(adjacencylist, f, w, "wsmaller_s", timelimit)

def dominated_model(adjacencylist, f, w, timelimit):
    return PCI_model_2(adjacencylist, f, w, "dominated", timelimit)

def dominated_model_weighted(adjacencylist, f, w, timelimit):
    return PCI_model_2(adjacencylist, f, w, "wdominated", timelimit)
