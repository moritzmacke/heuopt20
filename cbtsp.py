 
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
            raise ValueError("Invalid edge used in solution")
        
        if len(self.x) != self.inst.n:
            raise ValueError("Invalid length of solution")
        super().check()

    def construct(self, par, _result):
        """Scheduler method that constructs a new solution.

        #implment construction heuristic

        """
        if par == 0:
            self.initialize(par) # random order of nodes
        else:
            self.insert_const_heuristic(True)

    def insert_const_heuristic(self, select_random=False):
        """Starts with two nodes cycle and adds new nodes incrementally between
        those two that yields best objective function
        Does not work great the way it is, because a selected next point can
        force you to choose an invalid edge to insert in current tour...
        """
        
        k = 2
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
