"""A generic solution class for solutions that are represented by permutations of integers.
Adapted from pymhlib.PermutationSolution
"""

import numpy as np
from abc import ABC
from typing import List, Any, Tuple

from pymhlib.solution import VectorSolution, TObj
import random
from enum import IntEnum

class Step(IntEnum):
    RANDOM = 0
    BEST = 1
    FIRST = 2

class PermutationSolution(VectorSolution, ABC):
    """Solution that is represented by a permutation of 0,...length-1."""

    def __init__(self, length: int, init=True, **kwargs):
        """Initializes the solution with 0,...,length-1 if init is set."""
        super().__init__(length, init=False, **kwargs)
        if init:
            self.x[:] = np.arange(length)

    def copy_from(self, other: 'PermutationSolution'):
        super().copy_from(other)

    def initialize(self, k):
        """Random initialization."""
        np.random.shuffle(self.x)
        self.invalidate()

    def check(self):
        """Check if valid solution.

        :raises ValueError: if problem detected.
        """
        super().check()
        if set(self.x) != set(range(len(self.x))):
            raise ValueError("Solution is no permutation of 0,...,length-1")



    def neighborhood_search(self, generate_neighbors, apply_move_func, delta_function = None, strategy = Step.BEST) -> bool:
        """Generalized neighborhood search, completely seaches neighborhood generated by move generation function

        :param generate_neighbors: function used to generate the neighborhood, this is expected to explore the neighborhood
            in a random way, also gives inverse for each move
        :param apply_move_func: apply a move to solution (don't update or invalidate)
        :param delta_function: compute delta evaluation for move, if not given move is performed and reversered and 
            objective value is recomputed
        :param strategy:  Search until best solution or first improvement.

        :return: True if an improved solution has been found
        """
        
#        print("search with", generate_neighbors, apply_move_func, delta_function, strategy)

        #tried = set()

        n = self.inst.n
        best_delta = 0
        best_move = None

        for move, inverse in generate_neighbors(self):
            
#            print(move, inverse)
#            tried.add(move)
            if delta_function:
#                print("move", move, "inverse", inverse)
                
                delta = delta_function(self, *move)
                
#                delta2 = self.default_delta_eval(apply_move_func, move, inverse)
#                if delta != delta2:
#                    print("deltas",delta, delta2)
#                    assert(False)
                
            else:
                delta = self.default_delta_eval(apply_move_func, move, inverse)

            obj_new = self.obj_val + delta
            obj_best = self.obj_val + best_delta
            if self.is_better_obj(obj_new, obj_best) or strategy == Step.RANDOM:
                best_delta = delta
                best_move = move
                if strategy != Step.BEST:
                    break
        if best_move != None:
            apply_move_func(self, *best_move)
            self.obj_val += best_delta           
           
#            print(sorted(tried))
            return True
#        print(sorted(tried))
        return False

    """Generation methods for different neighborhood structures"""
    
    def generate_two_half_opt_neighborhood(self):
        """Generates two additional moves compared to 2opt
        moves are of the form (p,q,r): take rev sequence p..q (inclusive) and insert before r
        Whole 2.5opt functions currently quite slow...
        """
        n = self.inst.n
        order = np.arange(n)
        np.random.shuffle(order)
        for i,p1 in enumerate(order[:n-1]):
            for p2 in order[i+1:n]:
                p,q = (p1, p2) if (p1 < p2) else (p2, p1)
                #2opt move
                yield (p,q,q+1), (p,q,q+1)
                #other two
                yield (p,p,q+1), (q,q,p)
                yield (q,q,p), (p,p,q+1)
    
    
    
    def generate_two_exchange_neighborhood(self):
        n = self.inst.n
        order = np.arange(n)
        np.random.shuffle(order)
        for i,p1 in enumerate(order[:n-1]):
            for p2 in order[i+1:n]:
                p,q = (p1, p2) if (p1 < p2) else (p2, p1)
                yield (p,q), (p,q)
     
    def generate_two_opt_neighborhood(self):
        n = self.inst.n
        order = np.arange(n)
        np.random.shuffle(order)
        for i,p1 in enumerate(order[:n-1]):
            for p2 in order[i+1:n]:
                p,q = (p1, p2) if (p1 < p2) else (p2, p1)
                yield (p,q), (p,q)
                
    def generate_single_move_neighborhood(self):
        n = self.inst.n
        order = np.arange(n)
        np.random.shuffle(order)
        for i,p1 in enumerate(order[:n-1]):
            for p2 in order[i+1:n]:
                p,q = (p1, p2) if (p1 < p2) else (p2, p1)
                yield (p,q), (q+1,p) 
                               
                
    def generate_three_opt_neighborhood(self):
        raise NotImplementedError
    
    def generate_short_block_neighborhood(self, block_len = 3):
        """ """
        n = self.inst.n
        order = np.arange(n)
        np.random.shuffle(order)
        for i,p1 in enumerate(order[:n-1]):
            for p2 in order[i+1:n]:
                ins,blk = (p1, p2) if (p1 < p2) else (p2, p1)
                seqend = (blk + block_len) % n
                if seqend > blk or ins > seqend: 
                    #move seq at blk before ins
                    yield (ins,blk), (blk,ins) #cheating here and in apply function for reverse move...
                                               #just assume normal move is always p1<p2
#        raise NotImplementedError


    """Delta evaluation methods for different neighborhood structures"""
    
    def default_delta_eval(self, apply_move_func, move, inverse):
        """The solution is not changed.
        If not real delta function for a neighborhood is implemented.
        Here we perform the move, calculate the objective value from scratch and revert the move.
        """
        orig_obj = self.obj_val
        
#        print("before move", self.x)
        apply_move_func(self, *move)
#        print("after move ",self.x)
        self.invalidate()
        delta = self.obj() - orig_obj
        apply_move_func(self, *inverse)
#        print("after inv  ",self.x)
        self.obj_val = orig_obj
        return delta
    
    def single_move_delta_eval(self, p1: int, p2: int) -> int:
        """p2 moved in before p1"""
        
        n = len(self.x)
        d = self.inst.weights
        
        ipred_p1 = (p1 - 1) % n
        ipred_p2 = (p2 - 1) % n
        isucc_p2 = (p2 + 1) % n
        xpred_p1 = self.x[ipred_p1]
        xpred_p2 = self.x[ipred_p2]
        xsucc_p2 = self.x[isucc_p2]
        x_p1 = self.x[p1]
        x_p2 = self.x[p2]
        
        delta = d[xpred_p1][x_p2] + d[x_p2][x_p1] + d[xpred_p2][xsucc_p2] \
              - d[xpred_p1][x_p1] - d[xpred_p2][x_p2] - d[x_p2][xsucc_p2]
          
        return delta        
    
    def two_half_opt_move_delta_eval(self, p1: int, p2: int, p3: int) -> int:
        """ """
               
        n = len(self.x)
        d = self.inst.weights
        
        if p1 == p2:
            #single node move
            ipred_p1 = (p1 - 1) % n
            isucc_p1 = (p1 + 1) % n
            ipred_p3 = (p3 - 1) % n
            xpred_p1 = self.x[ipred_p1]
            xsucc_p1 = self.x[isucc_p1]
            xpred_p3 = self.x[ipred_p3]
            x_p1 = self.x[p1]
            x_p3 = self.x[p3 % n]
                
            if x_p1 == x_p3 or xsucc_p1 == x_p3:
                delta = 0
            else:
#                print("edges", (xpred_p1,xsucc_p1), (xpred_p3,x_p1), (x_p1,x_p3), "-", (xpred_p1,x_p1), (x_p1,xsucc_p1), (xpred_p3,x_p3))
                delta = d[xpred_p1][xsucc_p1] + d[xpred_p3][x_p1] + d[x_p1][x_p3] \
                        - d[xpred_p1][x_p1] - d[x_p1][xsucc_p1] - d[xpred_p3][x_p3]
            
        else:
            #normal 2opt move
            delta = self.two_opt_move_delta_eval(p1,p2)
            
        return delta      
        
    
    def two_exchange_move_delta_eval(self, p1: int, p2: int) -> int:
        """ This method performs the delta evaluation for exchanging p1 with p2.

        The function returns the difference in the objective function if the move would be performed,
        the solution, however, is not changed.
        """

        n = len(self.x)

        ipred_p1 = (p1 - 1) % n
        isucc_p1 = (p1 + 1) % n
        ipred_p2 = (p2 - 1) % n
        isucc_p2 = (p2 + 1) % n
        xpred_p1 = self.x[ipred_p1]
        xsucc_p1 = self.x[isucc_p1]
        xpred_p2 = self.x[ipred_p2]
        xsucc_p2 = self.x[isucc_p2]
        xp1 = self.x[p1]
        xp2 = self.x[p2]
        
        d = self.inst.weights
        delta = d[xpred_p1][xp2] + d[xp2][xsucc_p1] + d[xpred_p2][xp1] + d[xp1][xsucc_p2] \
                - d[xpred_p1][xp1] - d[xp1][xsucc_p1] - d[xpred_p2][xp2] - d[xp2][xsucc_p2]

        return delta    
    
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

    def three_opt_move_delta_eval(self, p1: int, p2: int, p3: int) -> int:
        """ Performs delta evaluation for moving the subsequence from p2 to p3 before p1 in self.x.

        The function returns the difference in the objective function if the move would be performed,
        the solution, however, is not changed.
        """
        assert (p1 < p2)
        assert (p2 < p3)
        n = len(self.x)

        x_p1_pred = self.x[(p1 - 1) % n]
        x_p1 = self.x[p1]
        x_p2_pred = self.x[(p2 - 1) % n]
        x_p2 = self.x[p2]
        x_p3 = self.x[p3]
        x_p3_succ = self.x[(p3 + 1) % n]
        d = self.inst_weights

        # Added edges: predecessor of p1 to p2, p3 to p1, predecessor of p2 to successor of p3
        # Lost edges: predecessor of p1 to p1, predecessor of p2 to p2, p3 to successor of p3
        delta = d[x_p1_pred][x_p2] + d[x_p3][x_p1] + d[x_p2_pred][x_p3_succ] \
            - d[x_p1_pred][x_p1] - d[x_p2_pred][x_p2] - d[x_p3][x_p3_succ]

        return delta

    def short_block_delta_eval(self, p1: int, p2: int) -> int:
        """ Performs delta evaluation for moving the subsequence of length 3 starting at p2 before p1 in self.x.

        The function returns the difference in the objective function if the move would be performed,
        the solution, however, is not changed.
        """
        assert (p1 < p2)
        n = len(self.x)

        x_p1_pred = self.x[(p1 - 1) % n]
        x_p1 = self.x[p1]
        x_p2_pred = self.x[(p2 - 1) % n]
        x_p2 = self.x[p2]
        x_subseq_end = self.x[(p2 + 2) % n]
        x_subseq_end_succ = self.x[(p2 + 3) % n]
        d = self.inst.weights

        # Added edges: predecessor of p1 to p2, subseq end to p1, predecessor of p2 to successor of subseq end
        # Lost edges: predecessor of p1 to p1, predecessor of p2 to p2, subseq end to successor of subseq end
        delta = d[x_p1_pred][x_p2] + d[x_subseq_end][x_p1] + d[x_p2_pred][x_subseq_end_succ] \
            - d[x_p1_pred][x_p1] - d[x_p2_pred][x_p2] - d[x_subseq_end][x_subseq_end_succ]

        return delta

    """Application methods for different neighborhood structures"""
    
    def apply_single_move(self, at: int, p: int):
        if at < p:
            self.x = np.concatenate((self.x[:at], self.x[p:p+1], self.x[at:p], self.x[p+1:]))
        else:
            self.x = np.concatenate((self.x[:p], self.x[p+1:at], self.x[p:p+1], self.x[at:]))
    
    def apply_two_half_opt_move(self, p1: int, p2: int, p3: int):
        
        if p1 == p2:
            if p3 > p1:
                self.x = np.concatenate((self.x[:p1], self.x[p1+1:p3], self.x[p1:p1+1], self.x[p3:]))
            elif p3 < p1:
                self.x = np.concatenate((self.x[:p3], self.x[p1:p1+1], self.x[p3:p1], self.x[p1+1:]))
            else:
                assert(False)
        else:
            self.apply_two_opt_move(p1,p2)

    def apply_two_exchange_move(self, p1: int, p2: int):
        """The values at positions p1 and p2 are exchanged in self.x.

        Note that the obj_val is not changed here nor is invalidate() called yet, as it is
        assumed that obj_val is updated by a corresponding delta evaluation.
        """
        x = self.x
        x[p1], x[p2] = x[p2], x[p1]
        
    def apply_two_opt_move(self, p1: int, p2: int):
        """The subsequence from p1 to p2 is inverted in self.x.

        Note that the obj_val is not changed here nor is invalidate() called yet, as it is
        assumed that obj_val is updated by a corresponding delta evaluation.
        """
        self.x[p1:(p2 + 1)] = self.x[p1:(p2 + 1)][::-1]

    def apply_three_opt_move(self, p1: int, p2: int, p3: int):
        """The subsequence from p2 to p3 is moved before p1 in self.x.

        Works the same way as apply_two_opt_move from the base class, so no value update or invalidation is done.
        """
        self.x = self.x[:p1] + self.x[p2:(p3+1)] + self.x[p1:p2] + self.x[(p3+1):]

    def _apply_short_block_move(self, p1: int, p2: int):
        """The subsequence of length 3 starting at p2 is moved before p1 in self.x.

        Works the same way as apply_two_opt_move from the base class, so no value update or invalidation is done.
        """
        self.x = self.x[:p1] + self.x[p2:(p2 + 3)] + self.x[p1:p2] + self.x[(p2 + 3):]
        


    def apply_short_block_move(self, p1: int, p2: int):
        """The subsequence of length 3 starting at p2 is moved before p1 in self.x.

        Works the same way as apply_two_opt_move from the base class, so no value update or invalidation is done.
        
        Only works for moves with p1 < p2, other case is assumed we want to reverse application of a move of first type
        """
        
        if p2 < p1:
            #special treatment for inverse move
       
            #here p1 is were original block was, p2 moved block location
            block_overflow = max((p1 + 3) - len(self.x), 0) 
            a_len = 3 - block_overflow
                        
            blk = p2-block_overflow
            block_a = self.x[blk:(blk+a_len)]
            block_b = self.x[(blk+a_len):blk+3]
             
            pre = self.x[:blk]
            mid = self.x[(blk + 3):(p1 + 3)]
            post = self.x[(p1 + 3):]
            
#            print(block_a, block_b, pre, mid, post)
            
            self.x = np.concatenate((block_b,pre,mid,block_a,post))
            
        else:
            #p1 < p2
            block_overflow = max((p2  + 3) - len(self.x), 0) 
            a_len = 3 - block_overflow
            block_a = self.x[p2:(p2 + a_len)]
            block_b = self.x[:block_overflow]
            
            pre = self.x[block_overflow:p1]
            mid = self.x[p1:p2]
            post = self.x[(p2 + 3):]
            
            self.x = np.concatenate((pre, block_a, block_b, mid, post))

        
    """Unchanged from pymhlib so far"""

    def random_two_exchange_move_delta_eval(self) -> Tuple[Tuple[int, int], TObj]:
        """Choose random move in the two-exchange neighborhood and perform delta eval., returning (p1, p2, delta_obj).

        The solution is not changed here yet.
        Primarily used in simulated annealing.
        """
        p1 = random.randint(0, len(self.x) - 2)
        p2 = random.randint(p1 + 1, len(self.x) - 1)
        delta_obj = self.two_exchange_move_delta_eval(p1, p2)
        return (p1, p2), delta_obj

    def random_two_opt_move_delta_eval(self) -> Tuple[Tuple[int, int], TObj]:
        """Choose random move in 2-opt neighborhood and perform delta evaluation, returning (move, delta_obj).

        The solution is not changed here yet.
        Primarily used in simulated annealing.
        """
        p1 = random.randrange(len(self.x)-1)
        p2 = random.randint(p1+1, len(self.x)-1)
        delta_obj = self.two_opt_move_delta_eval(p1, p2)
        return (p1, p2), delta_obj
    
    def partially_mapped_crossover(self, other: 'PermutationSolution') -> 'PermutationSolution':
        """Partially mapped crossover (PMX).

        Copies the current solution, selects a random subsequence from the other parent and realizes this subsequence
        in the child by corresponding pairwise exchanges.

        :param other: second parent
        :return: new offspring solution
        """

        size = len(self.x)

        # determine random subsequence
        begin = random.randrange(size)
        end = random.randrange(size - 1)
        if begin == end:
            end = end + 1
        if begin > end:
            begin, end = end, begin

        child = self.copy()

        # adopt subsequence from parent b by corresponding pairwise exchanges
        pos = np.empty(size, int)
        for i, elem in enumerate(child.x):
            pos[elem] = i
        for i in range(begin, end):
            elem = other.x[i]
            j = pos[elem]
            if i != j:
                elem_2 = child.x[i]
                child.x[i], child.x[j] = elem, elem_2
                pos[elem], pos[elem_2] = i, j
        child.invalidate()
        return child

    def cycle_crossover(self, other: 'PermutationSolution') -> 'PermutationSolution':
        """ Cycle crossover.

        A randomized crossover method that adopts absolute positions of the elements from the parents.

        :param other: second parent
        :return: new offspring solution
        """
        size = len(self.x)
        pos = np.empty(size, int)
        for i, elem in enumerate(self.x):
            pos[elem] = i

        # detect all cycles
        group = np.full(size, 0)
        group_id = 1
        for i in range(size):
            if group[i]:
                continue
            j = i
            while not group[j]:
                group[j] = group_id
                elem = other.x[j]
                j = pos[elem]
            group_id += 1

        # perform exchange
        child = self.copy()
        for i in range(size):
            if child.x[i] % 2 == 1:
                child.x[pos] = other.x[pos]
        child.invalidate()
        return child

    def edge_recombination(self, other: 'PermutationSolution') -> 'PermutationSolution':
        """ Edge recombination.

        This is a classical recombination operator for the traveling salesman problem, for example.
        It creates an adjacency list, i.e., a list of neighbors in the cyclically viewed parent permutations,
        for each element.
        A start element is randomly chosen.
        From this current element the next is iteratively determined by either choosing a neighbor with the smallest
        adjacency list (ties are broken randomly), or, if the list is of remaining neighbors is empty,
        by choosing some other not yet visited element at random.

        :param other: second parent
        :return new offspring solution
        """
        def append_if_not_contained(nbs, nb):
            if nb not in nbs:
                nbs.append(nb)
        size = len(self.x)
        adj_lists: List[List[int]] = [list() for _ in range(size)]
        for i, elem in enumerate(self.x):
            append_if_not_contained(adj_lists[elem], self.x[(i-1) % size])
            append_if_not_contained(adj_lists[elem], self.x[(i+1) % size])
        for i, elem in enumerate(other.x):
            append_if_not_contained(adj_lists[elem], other.x[(i-1) % size])
            append_if_not_contained(adj_lists[elem], other.x[(i+1) % size])
        unvisited = set(range(size))
        child = self.copy()
        elem = random.randrange(size)
        for i in range(size-1):
            # accept elem and remove it from unvisited and adjacency list
            child.x[i] = elem
            unvisited.remove(elem)
            for j in adj_lists[elem]:
                adj_lists[j].remove(elem)
            # select next elem
            if not adj_lists[elem]:
                sel = random.choice(list(unvisited))
            else:
                candidates = [adj_lists[elem][0]]
                degree = len(adj_lists[candidates[0]])
                for e2 in adj_lists[elem][1:]:
                    degree_e2 = len(adj_lists[e2])
                    if degree_e2 < degree:
                        candidates = [e2]
                    elif degree_e2 == degree:
                        candidates.append(e2)
                sel = random.choice(candidates)
                adj_lists[elem].clear()
            elem = sel
        child.x[-1] = elem
        child.invalidate()
        return child
    
    

