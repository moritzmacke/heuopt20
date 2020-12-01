from math import exp

import numpy as np
from pymhlib.sa import SA
from pymhlib.solution import TObj

class SA_CBTSP(SA):
    def metropolis_criterion(self, sol, delta_obj:TObj) -> bool:
        """Apply Metropolis criterion as acceptance decision determined by delta_obj and current temperature."""
        if sol.is_delta_improvement(delta_obj):
            return True
        return np.random.random_sample() <= exp(-abs(delta_obj) / self.temperature)