import random
import numpy as np
import math
from typing import Tuple, Any
from enum import IntEnum

import logging
from pymhlib.log import init_logger

from pymhlib.solution import TObj
from pymhlib.gvns import GVNS
from pymhlib.scheduler import Method
from pymhlib.solution import Solution

from permutation import PermutationSolution, Step
# from pymhlib.settings import get_settings_parser
from hamcycle import HamCycle
from greedyedge import GreedyEdgeConst


class Construct(IntEnum):
    NONE = 0
    GREEDY_EDGE = 1
    GREEDY_EDGE_RANDOM = 2
    HAMILTON_PATH = 3

class Neighbor(IntEnum):
    KOPT2 = 1
    KOPT3 = 2
    XCHG = 3        #swap two nodes
    SMOVE = 4       #move single node
    SBLOCK = 5      #move short sequence of nodes
    KOPT2HALF = 6

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
                edges.append((n1, n2, w))

        ws = sorted([e[2] for e in edges])
        minw = sum(ws[0:n])
        maxw = sum(ws[-n:])
        #M = maxw - min(minw,0)
        M = max(abs(minw),abs(maxw)) - sum(ws[0:n-1]) + 1
           
        # adjacency matrix
        weights = np.full((n, n), M)
        for (n1, n2, w) in edges:
#            print(n1,n2,w)
            weights[n1][n2] = weights[n2][n1] = w  #


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

    """Override comparision methods to deal with signed objective value"""

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


    def calc_objective(self):
        w = 0
        for i in range(self.inst.n - 1):
            w += self.inst.weights[self.x[i]][self.x[i + 1]]
        w += self.inst.weights[self.x[-1]][self.x[0]]
        return w  # abs(w)

    def check(self):
        """Check if valid solution.

        :raises ValueError: if problem detected.
        """

        # ?
        if self.obj() > self.inst.valid_threshold:
            invalid_edges = []
            for i in range(self.inst.n):
                j = (i + 1) % self.inst.n
                p, q = self.x[i], self.x[j]
                if self.inst.weights[p, q] == self.inst.bigM:
                    invalid_edges.append((p, q))
            raise ValueError("{} invalid edges used in solution: {}".format(len(invalid_edges), invalid_edges))

        if len(self.x) != self.inst.n:
            raise ValueError("Invalid length of solution")
        super().check()

    """Solution construction functions"""

    def construct(self, par, _result=None, alpha=0.1):
        """Scheduler method that constructs a new solution.
        """
        
        if par == Construct.GREEDY_EDGE:
            h = GreedyEdgeConst(self.inst.n, self.inst.edges)
            self.x[:] = h.construct(0)
            self.invalidate()
        elif par == Construct.GREEDY_EDGE_RANDOM:
            h = GreedyEdgeConst(self.inst.n, self.inst.edges)
            self.x[:] = h.construct(alpha)
            self.invalidate()
        elif par == Construct.HAMILTON_PATH:
            self.hamilton_const_heuristic()
        else:
            self.initialize(par)  # random order of nodes

    def hamilton_const_heuristic(self):
        """Try to find a hamiltonian cycle with valid edges, ignores edge weights
        """
        h = HamCycle(self.inst)
        self.x[:] = h.construct(True)
        self.invalidate()

    def shaking(self, par, result):
        """Scheduler method that performs shaking by 'par'-times swapping a pair of randomly chosen cities."""
        for _ in range(par):
            a = random.randint(0, self.inst.n - 1)
            b = random.randint(0, self.inst.n - 1)
            self.x[a], self.x[b] = self.x[b], self.x[a]
        self.invalidate()
        result.changed = True

    """Methods for local search"""
    
    def local_improve(self, _par, _result):
#        self.own_two_opt_neighborhood_search(True)

        neighbor, step = _par

        if neighbor == Neighbor.KOPT2:
#            self.own_two_opt_neighborhood_search(step == Step.BEST)
            gen = PermutationSolution.generate_two_opt_neighborhood
            app = PermutationSolution.apply_two_opt_move
            delta = PermutationSolution.two_opt_move_delta_eval
            self.neighborhood_search(gen, app, delta, step)
        elif neighbor == Neighbor.KOPT3:
            gen = PermutationSolution.generate_three_opt_neighborhood
            app = PermutationSolution.apply_three_opt_move
            delta = PermutationSolution.three_opt_move_delta_eval
            self.neighborhood_search(gen, app, delta, step)
        elif neighbor == Neighbor.XCHG:
            gen = PermutationSolution.generate_two_exchange_neighborhood
            app = PermutationSolution.apply_two_exchange_move
            delta = PermutationSolution.two_exchange_move_delta_eval
            self.neighborhood_search(gen, app, delta, step)
        elif neighbor == Neighbor.SMOVE:
            gen = PermutationSolution.generate_single_move_neighborhood
            app = PermutationSolution.apply_single_move
            delta = PermutationSolution.single_move_delta_eval
            self.neighborhood_search(gen, app, delta, step)
        elif neighbor == Neighbor.SBLOCK:
            gen = PermutationSolution.generate_short_block_neighborhood
            app = PermutationSolution.apply_short_block_move
            delta = PermutationSolution.short_block_delta_eval
            self.neighborhood_search(gen, app, delta, step)
        elif neighbor == Neighbor.KOPT2HALF:
            gen = PermutationSolution.generate_two_half_opt_neighborhood
            app = PermutationSolution.apply_two_half_opt_move
            delta = PermutationSolution.two_half_opt_move_delta_eval
            self.neighborhood_search(gen, app, delta, step)
        else:
            raise NotImplementedError

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
    
    def is_delta_improvement(self, delta):
        """Determines whether a given delta value is considered an improvement.
        A delta is an improvement if and only if adding it to the current objective value results in a value closer to 0.
        """
        return abs(self.obj_val + delta) < abs(self.obj_val)
    

if __name__ == '__main__':
    
    inst = CBTSPInstance("./instances/0010.txt")
    sol = CBTSPSolution(inst)
    ms = sorted([m for m in sol.generate_short_block_neighborhood()])
    print(ms)
    
    for (i,j),(ri,rj) in ms:
        sol2 = sol.copy()
        print(sol2)
        sol2.apply_short_block_move(i,j)
        print(sol2, (i,j))
        sol2.apply_short_block_move(ri,rj)
        print(sol2, (ri,rj),"\n")
        for i in range(len(sol.x)):
            assert(sol.x[i] == sol2.x[i])

