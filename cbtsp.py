 
import random
import numpy as np
import math
from typing import Tuple, Any

import logging
from pymhlib.log import init_logger

from pymhlib.permutation_solution import PermutationSolution
from pymhlib.solution import TObj
from pymhlib.gvns import GVNS
from pymhlib.scheduler import Method
from pymhlib.solution import Solution

#from pymhlib.settings import get_settings_parser
from hamcycle import HamCycle

class CBTSPInstance:
    """
    Attributes
        - n: 
        - weights: 
        - edges: 
        - valid_threshold
    """

    def __init__(self, file_name: str):
        """Read an instance from the specified file."""
        edges = []
        n = None
        m = None

        with open(file_name, "r") as f:
            lines = f.readlines()
            s = lines[0].split()
            n, m = int(s[0]), int(s[1])
            for line in lines[1:]:
                s = line.split()
                n1, n2, w = int(s[0]), int(s[1]), int(s[2])
                edges.append((n1,n2,w))
               
        ws = sorted([ e[2] for e in edges ])
        minw = sum(ws[0:n])
        maxw = sum(ws[-n:])
        M = maxw - min(minw,0)
            
        # adjacency matrix
        weights = np.full((n,n),M)
        for (n1, n2, w) in edges:
#            print(n1,n2,w)
            weights[n1][n2] = weights[n2][n1] = w #
            
        self.weights = weights
        self.n = n
        self.valid_threshold = maxw
        self.bigM = M
        self.edges = edges

    def __repr__(self):
        """Write out the instance data."""
        return f"n={self.n},\nweights={self.weights!r}\n"
    
    
class CBTSPSolution(PermutationSolution):
    """Solution to a TSP instance.

    Attributes
        - inst: associated TSPInstance
        - x: order in which cities are visited, i.e., a permutation of 0,...,n-1
    """

    to_maximize = False

    def __init__(self, inst: CBTSPInstance):
        super().__init__(inst.n, inst=inst)
        self.obj_val_valid = False

    def copy(self):
        sol = CBTSPSolution(self.inst)
        sol.copy_from(self)
        return sol

    def calc_objective(self):
        w = 0
        for i in range(self.inst.n - 1):
            w += self.inst.weights[self.x[i]][self.x[i + 1]]
        w += self.inst.weights[self.x[-1]][self.x[0]]
        return w #abs(w)

    def check(self):
        """Check if valid solution.

        :raises ValueError: if problem detected.
        """
        
        #?
        if self.obj() > self.inst.valid_threshold:
            invalid_edges = []
            for i in range(self.inst.n):
                j = (i+1) % self.inst.n
                p,q = self.x[i], self.x[j]
                if self.inst.weights[p,q] == self.inst.bigM:
                    invalid_edges.append((p,q))
            raise ValueError("{} invalid edges used in solution: {}".format(len(invalid_edges), invalid_edges))
        
        if len(self.x) != self.inst.n:
            raise ValueError("Invalid length of solution")
        super().check()

    def construct(self, par, _result):
        """Scheduler method that constructs a new solution.

        #implment construction heuristic

        """
        if par == 0:
            self.initialize(par) # random order of nodes
        elif par == 1:
            self.nn_const_heuristic_single(0)
        elif par == 2:
            self.edge_const_heuristic(False)
        elif par == 3:
            self.hamilton_const_heuristic()
        else:
            self.insert_const_heuristic(True)

    def hamilton_const_heuristic(self):
        """Try to find a hamiltonian cycle with valid edges, ignores edge weights
        """
        h = HamCycle(self.inst.n, self.inst.edges)
        self.x[:] = h.construct()
        self.invalidate()

    def insert_const_heuristic(self, select_random=False):
        """Starts with two nodes cycle and adds new nodes incrementally between
        those two that yields best objective function
        Does not work great the way it is, because a selected next point can
        force you to choose an invalid edge to insert in current tour...
        """
        
        w = self.inst.weights
        n = self.inst.n
        
        #start with shortest edge, is not necessesarily best though?
        edges = sorted(self.inst.edges, key = lambda e: abs(e[2]))
#        print(edges)
        p1 = edges[0][0]
        p2 = edges[0][1]
        newx = [p1,p2]
        ps = [p for p in range(n) if p not in newx]
        random.shuffle(ps)
        
#        print(ps)
        
        cur_obj = 2*w[p1][p2]
 #       print(cur_obj)
        
        for k in range(2,n):
            if select_random:
                p = ps.pop()
            else: #select farthest point for next insert
                #does not seem to work well though and is slow...
                farthest_w = 0
                farthest_p = None
                for p in ps:
                    best_w = min([ abs(w[p][q]) for q in newx])
                    if best_w > farthest_w:
                        farthest_w = best_w
                        farthest_p = p
#                    print("{} closest dist: {}".format(p, best_w))
                p = farthest_p
                ps.remove(p)
            best_j = 0
            p1 = newx[0]
            p2 = newx[1]
            best_delta = w[p1][p] + w[p][p2] - w[p1][p2]
#            print("obj insert {} between {} and {}: {}".format(p, p1, p2, best_delta+cur_obj))
            for j in range(1,k):
                p1 = newx[j]
                p2 = newx[(j+1)%k]
                delta = w[p1][p] + w[p][p2] - w[p1][p2]
#                print("obj insert {} between {} and {}: {}".format(p, p1,p2,delta+cur_obj))
                if abs(cur_obj+delta) < abs(cur_obj+best_delta):
                    best_j = j
                    best_delta = delta
            newx.insert(best_j+1, p)
            cur_obj += best_delta
#            print(newx)
        
        self.x[:] = newx
        self.invalidate()

    def nn_const_heuristic_single(self, starting_point):
        """Nearest neighbor-like heuristic, minimze cost of adding new point to tour.
        For single starting point, could repeat for every starting point
        """
        
        w = self.inst.weights
        n = self.inst.n
        
        x = [starting_point]
        ps = [q for q in range(n) if q != starting_point]
        
#        print(ps)
        
        cur_edge_sum = 0 #doesn't count return edge to close tour
 #       print(cur_obj)
        
        for k in range(1,n):
            p0, p = x[0], x[k-1]
            best_q = ps[0]
            best_new_obj = cur_edge_sum + w[p][best_q] + w[best_q][p0] 
            for q in ps:
                new_obj = cur_edge_sum + w[p][q] + w[q][p0]
#                print("add {} to tour: {}".format(q, new_obj))
                if abs(new_obj) < abs(best_new_obj):
                    best_q = q
                    best_new_obj = new_obj
            x += [best_q]
            cur_edge_sum += w[p][best_q]
            ps.remove(best_q)
#            print(x, cur_edge_sum + w[best_q][p0], ps)
        
        self.x[:] = x
        self.invalidate()

    def edge_const_heuristic(self, random_edge=False):
        """Greedily select edges without forming circles or trees (keep max degree <= 2)
        Seems like it produces fewer invalid edges in solution than nn for large instances?
        Could be used for randomized construction
        Could be sped up
        """
        n = self.inst.n
        edges = sorted(self.inst.edges, key=lambda e: e[2])
        setof = [i for i in range(n)]
        seledges = [[] for i in range(n)]
#        print(edges, setof, seledges)
        
        k = 0
        cur_obj = 0
        
        while edges:
            cand_edges = []
            for e in edges:
                p, q = e[0], e[1]
                if len(seledges[p]) > 1 or len(seledges[q]) > 1:
                    continue
                if setof[p] == setof[q] and k < n:
                    continue
                cand_edges.append(e)
#            print("filtered edges", cand_edges)
            if not cand_edges:
                break
            best_e = random.choice(cand_edges)
            best_obj = cur_obj + best_e[2]
            if not random_edge:
                for e in cand_edges:
                    new_obj = cur_obj + e[2]
                    if abs(new_obj) < abs(best_obj):
                        best_e = e
                        best_obj = new_obj
#            print("best edge", best_e)
            p,q = best_e[0], best_e[1]
            seledges[p].append(q)
            seledges[q].append(p)
            sq,sp = setof[q], setof[p]
            spq = min(sq,sp)
            for i in range(n):
                if setof[i] == sq or setof[i] == sp:
                    setof[i] = spq
#            print(["{}:{}".format(i,e) for (i,e) in enumerate(seledges)])
#            print(setof)
            edges = cand_edges
            cur_obj = best_obj
            k = k + 1
                
        
#        print(k)
        # construct tour from fragments, will have to add invalid edges
        
        # this whole procedure is pretty ugly and probably not correct...
        x = []
        visited = [False for i in range(n)]
        for i in range(n):
            if visited[i]:
                continue
            ps = [i]
            visited[i] = True
            if len(seledges[i]) >= 1:
                p = i
                q = seledges[p][0]
                ps.append(q)
                p = q
                while len(seledges[p]) > 1 and not visited[p]:
                    visited[p] = True
                    q1,q2 = seledges[p][0], seledges[p][1]
                    q = q1 if not visited[q1] else q2
                    if visited[q]:
                        break
                    ps.append(q)
                    p = q
                visited[p] = True
                q = seledges[p][0] #must be len(seledges[p] == 1 here, unless we have a full cycle?
                if not visited[q]:
                    ps.append(q)
                    visited[q] = True
                
                # search in other direction
                p = i  
                while len(seledges[p]) > 1:
                    q1,q2 = seledges[p][0], seledges[p][1]
                    q = q1 if not visited[q1] else q2
                    if visited[q]:
                        break
                    ps = [q] + ps
                    visited[q] = True
                    p = q
                q = seledges[p][0]
                if not visited[q]:
                    ps = [q] + ps
                    visited[q] = True
                    
#            print(ps)
            x += ps

        
#        print(x, [i for i in x if x.count(i) > 1])
        self.x[:] = x
        self.invalidate()

 

    def is_better(self, other: "Solution") -> bool:
        """Return True if the current solution is better in terms of the objective function than the other."""
        return abs(self.obj()) < abs(other.obj()) 

    def is_worse(self, other: "Solution") -> bool:
        """Return True if the current solution is worse in terms of the objective function than the other."""
        return abs(self.obj()) > abs(other.obj()) 

    @classmethod
    def is_better_obj(cls, obj1: TObj, obj2: TObj) -> bool:
        """Return True if the obj1 is a better objective value than obj2."""
        return abs(obj1) < abs(obj2)

    @classmethod
    def is_worse_obj(cls, obj1: TObj, obj2: TObj) -> bool:
        """Return True if obj1 is a worse objective value than obj2."""
        return abs(obj1) > abs(obj2)









#
    def shaking(self, par, result):
        """Scheduler method that performs shaking by 'par'-times swapping a pair of randomly chosen cities."""
        for _ in range(par):
            a = random.randint(0, self.inst.n - 1)
            b = random.randint(0, self.inst.n - 1)
            self.x[a], self.x[b] = self.x[b], self.x[a]
        self.invalidate()
        result.changed = True

    def local_improve(self, _par, _result):
        self.own_two_opt_neighborhood_search(True)
        
    #
    def own_two_opt_neighborhood_search(self, best_improvement) -> bool:
        """Systematic search of the 2-opt neighborhood, i.e., consider all inversions of subsequences.

        The neighborhood is searched in a randomized ordering.
        
        :param best_improvement:  if set, the neighborhood is completely searched and a best neighbor is kept;
            otherwise the search terminates in a first-improvement manner, i.e., keeping a first encountered
            better solution.

        :return: True if an improved solution has been found
        """
               
        n = self.inst.n
        best_delta = 0
        best_p1 = None
        best_p2 = None
        order = np.arange(n)
        np.random.shuffle(order)
        for idx, p1 in enumerate(order[:n - 1]):
            for p2 in order[idx + 1:]:
                if p1 > p2:
                    p1, p2 = p2, p1
                # consider the move that self.x from position p1 to position p2
                delta = self.two_opt_move_delta_eval(p1, p2)
                obj_new = self.obj_val + delta
                obj_best = self.obj_val + best_delta
                if self.is_better_obj(obj_new, obj_best):
                    if not best_improvement:
                        self.apply_two_opt_move(p1, p2)
                        self.obj_val += delta
                        return True
                    best_delta = delta
                    best_p1 = p1
                    best_p2 = p2
        if best_p1:
            self.apply_two_opt_move(best_p1, best_p2)
            self.obj_val += best_delta
            return True
        return False

    def two_opt_move_delta_eval(self, p1: int, p2: int) -> int:
        """ This method performs the delta evaluation for inverting self.x from position p1 to position p2.

        The function returns the difference in the objective function if the move would be performed,
        the solution, however, is not changed.
        """
        assert (p1 < p2)
        n = len(self.x)
        if p1 == 0 and p2 == n - 1:
           # reversing the whole solution has no effect
            return 0
        prev = (p1 - 1) % n
        nxt = (p2 + 1) % n
        x_p1 = self.x[p1]
        x_p2 = self.x[p2]
        x_prev = self.x[prev]
        x_next = self.x[nxt]
        d = self.inst.weights
        delta = d[x_prev][x_p2] + d[x_p1][x_next] - d[x_prev][x_p1] - d[x_p2][x_next]

        return delta

    def random_move_delta_eval(self) -> Tuple[Any, TObj]:
        """Choose a random move and perform delta evaluation for it, return (move, delta_obj)."""
        return self.random_two_opt_move_delta_eval()

    def apply_neighborhood_move(self, move):
        """This method applies a given neighborhood move accepted by SA,
            without updating the obj_val or invalidating, since obj_val is updated incrementally by the SA scheduler."""
        self.apply_two_opt_move(*move)

    def crossover(self, other: 'CBTSPSolution') -> 'CBTSPSolution':
        """Perform edge recombination."""
        return self.edge_recombination(other)
