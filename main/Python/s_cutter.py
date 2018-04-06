import cplex
import random

random.seed(314)

class S_cutter:

    # -----CONSTRUCTOR
    def __init__(self, adjacencylist, f):
        self.adjacencylist = adjacencylist
        self.f = f
        self.N = len(f)
        self.v_counter = [0] * len(f)
        self.variables_used = [0] * len(f)
        self.counter = 0


#############################################################################
    # -----AUXILIARY METHODS
    def infect_graph(self, infected):
        adjacencylist = self.adjacencylist
        f = self.f
        V = range(len(f))
        new_f = [_ for _ in f]
        neighboors_infected = [0 for v in V]
        queue = [v for v in V if infected[v]]
        while queue:
            v = queue.pop(0)
            new_f[v] = -2
            for u in adjacencylist[v]:
                if infected[u]:
                    continue
                neighboors_infected[u] += 1
                new_f[u] -= 1
                if new_f[u] <= 0:
                    queue.append(u)
                    infected[u] = True
        #print(new_f)
        return new_f

    # finds | N(v) intersect (V/S)|
    def neighbors_outside_s (self, S, v):
        solution = 0
        S_bool = [False for _ in range(self.N)]
        for s in S:
            S_bool[s] = True
        for u in self.adjacencylist[v]:
            if not S_bool[u]:
                solution += 1
        #print(S, v, solution)
        return solution

    # for a given S, tries to find violated inequalities
    def S_induced_inequality (self, S):
        # generate every possible S
        # there is a violated inequality if min (c0-1(S) - |N(v) intersection (V/S))
        assert(S)
        f = self.f
        min_S_infected = f[S[0]] + 1
        min_vertex = -1
        # finds min(f(v) - |N_G(v)) inter (V/S)|
        for v in S:
            new_min_val = f[v] - self.neighbors_outside_s(S, v)
            if (new_min_val < min_S_infected):
                min_S_infected = new_min_val
                min_vertex = v
        # a vertex should be choosen as min
        assert(min_vertex != -1)
        return min_S_infected

    # given a set of infected vertices and an f, creates constraints based on
    # the S set of remaining vertices
    def finds_S_constraints(self, infected_vertices, new_f):
        adjacencylist = self.adjacencylist
        f = self.f
        V = range(len(new_f))
        visited = [(new_f[v] <= 0) for v in V]
        components_limits = []
        # finds connected component
        for w in V:
            if visited[w]: continue
            queue = [w]
            visited[w] = True
            S = []
            while queue:
                v = queue.pop(0)
                if new_f[v] > 0:
                    S.append(v)
                    for u in adjacencylist[v]:
                        if not visited[u]:
                            queue.append(u)
                            visited[u] = True
            if S:
                S_lower_bound = self.S_induced_inequality(S)

                if S_lower_bound > 0:
                    #assignment_constraint = cplex.SparsePair (ind = S,
                    #     val = [1] * len(S))
                    components_limits.append([S, S_lower_bound])

        self.constraints = [_ for _ in components_limits]
        return (len(components_limits) != 0)


############################################################################
    # -------- WEIGHTED SHUFFLE:

    # Select random vertex based on how many variables were used
    # the vector counting the use sums to 1
    def recursive_finds_rw_vertex(self, i, j, r, v_counter):
        assert(i < j) # asserts the range is [i, .. , j)
        mid = (i + j) / 2
        # if mid is outside the vector range, found the solution
        if mid - 1 < i or mid + 1 == j:
            return mid
        assert (i <= mid and mid < j)
        if v_counter[mid] <= r:
            if mid + 1 == j  or r < v_counter[mid + 1]:
                return mid
            else:
                return self.recursive_finds_rw_vertex(mid, j, r, v_counter)

        return self.recursive_finds_rw_vertex(i, mid, r, v_counter)

    # Finds an order to select vertices in a weighted but randomized matter
    def weighted_shuffle(self):
        weights_vector = [_ for _ in self.v_counter]
        shuffle_positions = range(self.N)
        for i in range(self.N):
            r = random.uniform(0, 1)
            #shuffles from position i to N, because 0 to i-1 is shuffled
            j = self.recursive_finds_rw_vertex(i, self.N, r, weights_vector)
            #found he next value in j to the shuffled vector, puts in position i
            weights_vector[i], weights_vector[j] = weights_vector[j], weights_vector[i]
            shuffle_positions[i], shuffle_positions[j] = shuffle_positions[j], shuffle_positions[i]
        return shuffle_positions

    def shuffle_positions_in_range(self, i, j, weights_vector):
        positions = [weights_vector[p][0] for p in range(i, j)]
        random.shuffle(positions)
        for p in range(i, j):
            weights_vector[p][0] = positions[p - i]

    def shuffle_equal_positions(self, weights_vector):
        current_weight = weights_vector[0][1]
        i = 0
        for j in range(1, self.N):
            if current_weight != weights_vector[j][1]:
                self.shuffle_positions_in_range(i, j, weights_vector)
                i = j
                current_weight = weights_vector[j][1]
        self.shuffle_positions_in_range(i, self.N, weights_vector)
        #asserts everything is in order
        for i in range(1, self.N):
            assert(weights_vector[i - 1][1] >= weights_vector[i][1])

    def weighted_positions(self):
        weights_vector = [[i, self.v_counter[i]] for i in range(self.N)]
        weights_vector = sorted(weights_vector, key = lambda weight: weight[1], reverse=True)
        #print(weights_vector)
        self.shuffle_equal_positions(weights_vector)
        #print("NEWWWWWWW ", weights_vector)
        return [t[0] for t in weights_vector]

    def normal_shuffle(self):
        r = range(self.N)
        random.shuffle(r)
        return r

    def select_next_vertex(self, vertices, infected, position):
        for i in range(position, len(vertices)):
            if not infected[vertices[i]]:
                position = i
                return vertices[i]
        return -1

    # When new constraints are found, reweight the vector that counts variables
    # used in the constraints
    def reweight_vector_vertices_selected(self, bounds):
        for constraint in bounds:
            for i in constraint[0]:
                self.variables_used[i] += 1
                self.counter += 1
        # if counters gets to big, reduces it (but we will have bigger problems)
        if self.counter > 100000:
            for i in range(self.N):
                self.variables_used[i] = self.variables_used[i] / 1000
        # updates v_counter
        if self.counter == 0:
            return
        for i in range(self.N):
            self.v_counter[i] = float(self.variables_used[i]) / self.counter
        #print("HEEEERE !!!!constraint ", self.variables_used)

    def weighted_option_selected(self, option):
        if option[0] == 'w':
            #return self.weighted_shuffle()
            return self.weighted_positions()
        return self.normal_shuffle()

#############################################################################
    #------ FIRST SIMPLE S MODEL

    # gets the not infected vertices, uses that as s and creates constraints
    def finds_s_model_constraints(self, infected_vertices, option):
        self.new_f = self.infect_graph(infected_vertices)
        self.infected_vertices = infected_vertices
        return self.finds_S_constraints(infected_vertices, self.new_f)

##############################################################################
    # ------ SMALLER S
    # NEW CUTTING ALGORITHM THAT INFECTS THE REST and gets smaller S

    # given a range and a set of forbidden vertices, selects the one avaliable
    # if no vertex is selected, returns -1
    def select_random_vertex(self, infected):
        num_not_infected = 0
        for i in infected:
            if (not i): num_not_infected += 1
        r = random.randint(0, num_not_infected - 1)
        j = 0
        for i in range(len(infected)):
            if (infected[i]):
                continue
            if (j == r):
                return i
            j += 1
        return -1

    # given one new infected vertex in the graph, carries on the infection
    # this vertex cannot be already infected
    def infect_one_vertex(self, vi, infected, new_f):
        adjacencylist = self.adjacencylist
        queue = [vi]
        assert(not infected[vi])
        infected[vi] = True
        new_f[vi] = -2
        while queue:
            v = queue.pop(0)
            for u in adjacencylist[v]:
                if not infected[u]:
                    new_f[u] -= 1
                    if new_f[u] <= 0:
                        queue.append(u)
                        infected[u] = True
                        new_f[u] = -2


    # modifies infected vertices, so everything outside one component counts as
    # infected. Return if we still have components
    def can_select_random_component(self, infected):
        cv = -1
        for v in range(len(infected)):
            if not infected[v]:
                cv = v
        if cv == -1:
            return False
        queue = [cv]
        component = [False for _ in infected]
        while queue:
            v = queue.pop(0)
            component[v] = True
            for u in self.adjacencylist[v]:
                if not infected[u] and not component[u]:
                    component[u] = True
                    queue.append(u)
        # modifies infected vertices, so if v is outside component, treat as
        # infect
        infected = [(infected[v] or not component[v]) for v in range(len(infected))]
        return True


    #this function returns if we found constraints for a smaller S
    def finds_s_smaller_constraints(self, infected, option):
        already_infected = [_ for _ in infected]
        new_f = self.infect_graph(already_infected)
        self.finds_S_constraints(already_infected, new_f)
        bounds = self.constraints
        # gets a vector with the section order of vertices based on the option
        # weighted or not
        selection_order = self.weighted_option_selected(option)
        #print(selection_order)
        position = 0
        #keeps selecting one of the avaliable components
        while (self.can_select_random_component(already_infected)):
            #v = self.select_random_vertex(already_infected)
            v = self.select_next_vertex(selection_order, already_infected, position)
            assert(v != -1)
            self.infect_one_vertex(v, already_infected, new_f)
            if (self.finds_S_constraints(already_infected, new_f)):
                bounds = self.constraints
        if option[0] == 'w':
            self.reweight_vector_vertices_selected(bounds)
            #print("v_vector ", self.variables_used )

        #print("Final bounds ", bounds)
        self.constraints = bounds
        return len(bounds) != 0

##############################################################################
    # --------- DOMINATION
    def constraint_contained_in_other(self, constraint, bigger_constraint):
        #print("constraint ", constraint[0])
        assert(len(constraint[0]) == self.N)
        for i in range(self.N):
            if constraint[0][i] == 1 and bigger_constraint[0][i] == 0:
                return False
        return True

    def constraint_is_dominated(self, c, bc):
        constraint = [[0] * self.N, c[1]]
        bigger_constraint = [[0] * self.N, bc[1]]
        for i in range(len(c[0])):
            constraint[0][c[0][i]] = 1
        for i in range(len(bc[0])):
            bigger_constraint[0][bc[0][i]] = 1
        if not self.constraint_contained_in_other(constraint, bigger_constraint):
            return False
        left_difference = 0
        for i in range(self.N):
            if bigger_constraint[0][i] == 1:
                # number of different variables between each constraint
                left_difference += 1 - constraint[0][i]
        # difference in right side of constraints
        right_difference = bigger_constraint[1] - constraint[1]
        # if there are to many variables removed, compared with the min number
        # of activated vertices mandatory, then it is not dominated and must be
        # added this new constraint
        return left_difference <= right_difference

    def add_not_dominated_new_constraints(self, bounds):
        new_add_bounds = []
        for constraint in self.constraints:
            constraint_not_dominated = True
            for bigger_constraint in bounds:
                if self.constraint_is_dominated(constraint, bigger_constraint):
                    constraint_not_dominated = False
                    break
            if constraint_not_dominated:
                new_add_bounds.append(constraint)
        bounds += new_add_bounds

    def add_not_dominated_new_constraints_smaller(self, bounds):
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


    #this function returns if we found constraints for a smaller S and domination
    def finds_s_with_domination_constraints(self, infected, option):
        already_infected = [_ for _ in infected]
        new_f = self.infect_graph(already_infected)
        self.finds_S_constraints(already_infected, new_f)
        bounds = self.constraints
        # gets a vector with the section order of vertices based on the option
        # weighted or not
        selection_order = self.weighted_option_selected(option)
        #print(selection_order)
        position = 0
        #keeps selecting one of the avaliable components
        while (self.can_select_random_component(already_infected)):
            #v = self.select_random_vertex(already_infected)
            v = self.select_next_vertex(selection_order, already_infected, position)
            assert(v != -1)
            self.infect_one_vertex(v, already_infected, new_f)
            if (self.finds_S_constraints(already_infected, new_f)):
                self.add_not_dominated_new_constraints(bounds)

        if option[0] == 'w':
            self.reweight_vector_vertices_selected(bounds)
            #print("v_vector ", self.variables_used )

        #print("Final bounds ", bounds)
        self.constraints = bounds
        return len(bounds) != 0

#############################################################################

    def finds_constraints(self, infected, option):
        copy_infected = [_ for _ in infected]
        cases = { "s_model": self.finds_s_model_constraints,
            "smaller_s": self.finds_s_smaller_constraints,
            "wsmaller_s": self.finds_s_smaller_constraints,
            "dominated": self.finds_s_with_domination_constraints,
            "wdominated": self.finds_s_with_domination_constraints
        }
        #print(cases[option])
        return cases[option](copy_infected, option)
