import random
import numpy as np
import math
import sys
from typing import Tuple, Any
from enum import IntEnum

import logging
from pymhlib.log import init_logger

from pymhlib.solution import TObj
from pymhlib.gvns import GVNS
from pymhlib.scheduler import Method
from pymhlib.solution import Solution

from permutation import PermutationSolution, Step, NeighborhoodSpec
# from pymhlib.settings import get_settings_parser
from hamcycle import HamCycle
from greedyedge import GreedyEdgeConst


class Construct(IntEnum):
    NONE = 0
    GREEDY_EDGE = 1
    GREEDY_EDGE_RANDOM = 2
    HAMILTON_PATH = 3

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
        self.valid_threshold = max(abs(minw), abs(maxw))
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

    def get_invalid_edges(self):
        """Get the list of invalid edges in solution.
        """
        invalid_edges = []
        for i in range(self.inst.n):
            j = (i + 1) % self.inst.n
            p, q = self.x[i], self.x[j]
            if self.inst.weights[p, q] == self.inst.bigM:
                invalid_edges.append((p, q))

        return invalid_edges

    def check(self):
        """Check if valid solution.

        :raises ValueError: if problem detected.
        """
        invalid_edges = self.get_invalid_edges()
        if len(invalid_edges) > 0:
            print("warning: {} invalid edges used in solution: {}".format(len(invalid_edges), invalid_edges), file=sys.stderr)

        if len(self.x) != self.inst.n:
            raise ValueError("Invalid length of solution")
        super().check()

    """Solution construction functions"""

    def construct(self, par, _result=None, alpha=0.1, timeout=20):
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
            """Try to find a hamiltonian cycle with valid edges, ignores edge weights"""
            h = HamCycle(self.inst)
            self.x[:] = h.construct(True, timeout)
            self.invalidate()
        else:
            self.initialize(par)  # random order of nodes


    def random_construct(self, par, _result=None):
               
        alpha_step = 1.0/len(self.inst.edges)
        alpha_val = max(par['alpha'], alpha_step)
        min_alpha = 2*alpha_step       

        self.construct(Construct.GREEDY_EDGE_RANDOM, alpha=alpha_val)
#        print(alpha_val, self.obj())
                   
        alpha_val += max(alpha_step, 0.01)
        if alpha_val > min(1.0, (alpha_step* max(100,len(self.inst.edges)/10) )):
            alpha_val = min_alpha
            
        par['alpha'] = alpha_val
        

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
        neighbor, step = _par
        if isinstance(neighbor, NeighborhoodSpec):
            self.neighborhood_search(neighbor, step)
        else:
            raise NotImplementedError
    
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
    
