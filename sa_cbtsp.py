from math import exp
from typing import List, Callable

import numpy as np
from pymhlib.sa import SA
from pymhlib.scheduler import Method
from pymhlib.solution import TObj, Solution


class SA_CBTSP(SA):
    """
    Adapted Simulated Annealing metaheuristic for CBTSP.

    Added attributes:
    -) reheat_iter: Number of iterations since current best was found that have to pass to trigger a reheat
    -) reheat_temp: Temperature that is reheated to
    """

    def __init__(self, sol: Solution, meths_ch: List[Method], random_move_delta_eval: Callable,
    apply_neighborhood_move: Callable, iter_cb: Callable, own_settings: dict = None, consider_initial_sol = False):
        """Initialization.

        :param sol: solution to be improved
        :param meths_ch: list of construction heuristic methods
        :param random_move_delta_eval: function that chooses a random move and determines the delta in the obj_val
        :param apply_neighborhood_move: apply neighborhood move method return by propose method
        :param iter_cb: callback for each iteration passing iteration number, proposed sol, accepted sol, temperature,
            and acceptance
        :param own_settings: optional dictionary with specific settings
        :param consider_initial_sol: if true consider sol as valid solution that should be improved upon; otherwise
            sol is considered just a possibly uninitialized of invalid solution template
        """
        super().__init__(sol, meths_ch, random_move_delta_eval, apply_neighborhood_move, iter_cb, own_settings, consider_initial_sol)
        self.reheat_iter = self.own_settings.mh_sa_reheat_iter
        self.reheat_temp = self.own_settings.mh_sa_T_reheat
        self.last_reheat_iteration = 0

    def metropolis_criterion(self, sol, delta_obj:TObj) -> bool:
        """Apply Metropolis criterion as acceptance decision determined by delta_obj and current temperature."""
        if sol.is_delta_improvement(delta_obj):
            return True
        return np.random.random_sample() <= exp(-abs(delta_obj) / self.temperature)

    def cool_down(self):
        """Apply geometric cooling."""
        self.temperature *= self.own_settings.mh_sa_alpha
        # Reheat if necessary
        if self.iteration - max(self.incumbent_iteration, self.last_reheat_iteration) > self.reheat_iter:
            self.temperature = self.reheat_temp
            self.last_reheat_iteration = self.iteration