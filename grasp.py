"""Simplified from pymhlib.GVNS 
"""

from typing import List
import time

from pymhlib.scheduler import Method, Scheduler, Result
from pymhlib.settings import get_settings_parser
from pymhlib.solution import Solution


parser = get_settings_parser()


class GRASP(Scheduler):
    """.

    Attributes
        - sol: solution object, in which final result will be returned
        - meth_ch: construction heuristic method
        - meth_li: local improvement method
    """

    def __init__(self, sol: Solution, meth_ch: Method, meth_li: Method, own_settings: dict = None, consider_initial_sol=False):
        """Initialization.

        :param sol: solution to be improved
        :param meth_ch: construction heuristic method
        :param meth_li: local improvement method
        :param own_settings: optional dictionary with specific settings
        :param consider_initial_sol: if true consider sol as valid solution that should be improved upon; otherwise
            sol is considered just a possibly uninitialized of invalid solution template
        """
        super().__init__(sol, [meth_ch,meth_li], own_settings, consider_initial_sol)
        self.meth_ch = meth_ch
        self.meth_li = meth_li

    #just search till local optimum? 
    def localsearch(self, sol: Solution) -> Result:
        """local search until no better solution found or termination condition

        :returns: 
        """
        sol2 = sol.copy()
        while True:
#            print("improve", sol2, sol2.obj())
            res = self.perform_method(self.meth_li, sol2)
            if sol2.is_better(sol):
#                print("improved", sol2.obj())
                sol.copy_from(sol2)
                if res.terminate:
                    return res
            else:
#                print("couldn't be improved", sol2, sol2.obj())
                return res


    def grasp(self):
        """Perform GRASP to given solution."""
        while True:
            sol = self.incumbent.copy()
            res = self.perform_method(self.meth_ch, sol)
#            print("constructed", sol.obj())
            if res.terminate:
                break
            
            res = self.localsearch(sol)
#            print("localsearch", sol.obj())
            if sol.is_better(self.incumbent):
                self.update_incumbent(sol, time.process_time() - self.time_start)
            if res.terminate:
                break
            


    def run(self) -> None:
        """ """
        assert self.incumbent_valid or self.meth_ch
        self.grasp()
