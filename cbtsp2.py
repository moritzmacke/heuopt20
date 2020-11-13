 
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


        Here we just call initialize.
        """
        self.initialize(par)




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

if __name__ == '__main__':
    from pymhlib.settings import parse_settings
    parse_settings()
    
    inst = CBTSPInstance("instances/0010.txt")
    
        
    print(inst)
    
    sol = CBTSPSolution(inst)   
    print(sol.obj())
    
    
    
    init_logger()
    logger = logging.getLogger("pymhlib")
    logger.info("...")
    
    
    
#    (self, sol: Solution, meths_ch: List[Method], meths_li: List[Method], meths_sh: List[Method],
#                 own_settings: dict = None, consider_initial_sol=False):
    alg = GVNS(sol, [Method(f"ch0", CBTSPSolution.construct, 0)], [Method(f"li1", CBTSPSolution.local_improve, 1)], [Method(f"sh5", CBTSPSolution.shaking, 5)], {'mh_checkit': False, 'mh_titer': -100, 'mh_tciter': -1, 'mh_ttime': 10, 'mh_tctime': -1, 'mh_tobj': -1, 'mh_lnewinc': True, 'mh_lfreq': 0} )
    
    alg.run()
    
    
    print(sol)
    print(sol.obj())
    
    alg.method_statistics()
    alg.main_results()
    
#    run_optimization('CBTSP', CBTSPInstance, CBTSPSolution, "instances/0010.txt")
