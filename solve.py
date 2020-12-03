import logging
from typing import Callable
import time

from pymhlib.gvns import GVNS
from pymhlib.log import init_logger
from pymhlib.scheduler import Method
from pymhlib.settings import parse_settings, settings, get_settings_parser, get_settings_as_str
from pymhlib.solution import Solution

from cbtsp import *
from grasp import GRASP
from sa_cbtsp import SA_CBTSP

inst_dir = "instances/"

def add_general_arguments_and_parse_settings(default_inst_file: str = '0010.txt', seed: int = 0):
    """Some general parameters are registered and the settings are parsed.

    :param seed: optional seed value for the random number generators; 0: random initialization
    :param default_inst_file: default instance file to be loaded and solved
    """
    parser = get_settings_parser()
    parser.add_argument("--alg", type=str, default='vnd', help='optimization algorithm to be used '
                                                                '(just_const, just_rconst, grasp, lsearch, gvns)')
    parser.add_argument("--inst_file", type=str, default=default_inst_file,
                        help='problem instance file')
    parser.add_argument("--use_delta_eval", type=bool, default=True, help="")
    parser.add_argument("--neighborhood", type=str, default="2opt", help="(2opt, ...)")
    parser.add_argument("--step_function", type=str, default="best", help="(best, first, random) improvement")
    parse_settings(seed=seed)


if __name__ == '__main__':
    
    add_general_arguments_and_parse_settings('0010.txt', 0)

    init_logger()
    logger = logging.getLogger("pymhlib")
    logger.info(get_settings_as_str())
    
    ownsettings = {'mh_titer': -1, 'mh_workers': 1, 'mh_ttime': 20 if settings.mh_ttime < 0 else settings.mh_ttime}
    logger.info(ownsettings)
    
    
    instance = CBTSPInstance(inst_dir + settings.inst_file)
    logger.info("instance read:\n" + str(instance))
    solution = CBTSPSolution(instance)
    
    
    
    if settings.alg == 'just_const':
        solution.construct(Construct.HAMILTON_PATH, None)
        solution.check()
        print("obj", solution.obj())
    elif settings.alg == 'just_rconst':
        best_sol = solution.copy()
        start = time.process_time()
        elapsed = 0
        found_obj_vals = set() # this will grow too large eventually, need better way to see if should increase alpha
        alpha_val = 0
        while elapsed < ownsettings['mh_ttime']:
            solution.construct(Construct.GREEDY_EDGE_RANDOM, alpha=alpha_val)
            elapsed = time.process_time() - start
            print(elapsed, "alpha", alpha_val, "obj", solution.obj(), "best", best_sol.obj())
#            print(len(found_obj_vals))
            if solution.obj() in found_obj_vals:
                alpha_val = min(alpha_val + 0.01, 0.5)
            found_obj_vals.add(solution.obj())
            if solution.is_better(best_sol):
                best_sol, solution = solution, best_sol

        print("best obj", best_sol.obj())
        best_sol.check()

    elif settings.alg == 'grasp':
        alg = GRASP(solution, Method("rconst", CBTSPSolution.construct, Construct.GREEDY_EDGE_RANDOM), Method("search", CBTSPSolution.local_improve, (Neighbor.KOPT2, Step.BEST)), ownsettings)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
    elif settings.alg == 'lsearch':
        
#        print(solution, solution.obj())
#        solution.local_improve((Neighbor.KOPT2, Step.BEST),None)
#        print(solution, solution.obj())
        
        raise NotImplementedError
    elif settings.alg == 'gvns':
        alg = GVNS(solution, [Method(f"ch0", CBTSPSolution.construct, Construct.HAMILTON_PATH)], [Method("li_2opt_best", CBTSPSolution.local_improve, (Neighbor.KOPT2, Step.BEST))], [Method(f"sh{i}", CBTSPSolution.shaking, i) for i in range(1, 2)], ownsettings)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
    elif settings.alg == "sa":
        sa_settings = {
            'mh_titer': -1,
            'mh_tciter': 10000, # Shortcut: Abort after 10000 non-improving iterations - remove for real tests
            'mh_ttime': 15*60 # Limited to 15 min CPU time
        }
        alg = SA_CBTSP(solution, [Method("rconst", CBTSPSolution.construct, Construct.GREEDY_EDGE_RANDOM)], CBTSPSolution.random_move_delta_eval, CBTSPSolution.apply_neighborhood_move, None, sa_settings)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
    elif settings.alg == "vnd":
        # Which neighborhoods do we want to use here?
        alg = GVNS(solution, [], random.sample([
            Method("li_2opt_best", CBTSPSolution.local_improve, (Neighbor.KOPT2, Step.BEST)),
            Method("li_3opt_best", CBTSPSolution.local_improve, (Neighbor.KOPT3, Step.BEST)),
            Method("li_xchg_best", CBTSPSolution.local_improve, (Neighbor.XCHG, Step.BEST)),
            Method("li_smove_best", CBTSPSolution.local_improve, (Neighbor.SMOVE, Step.BEST)),
            Method("li_sblock_best", CBTSPSolution.local_improve, (Neighbor.SBLOCK, Step.BEST)),
            Method("li_2.5opt_best", CBTSPSolution.local_improve, (Neighbor.KOPT2HALF, Step.BEST))
        ], 6), [], ownsettings)
        alg.vnd(solution.copy()) #Needs to be a copy or update_incumbent in algorithm will override the argument
        logger.info("")
        alg.method_statistics()
        alg.main_results()

#    print(solution, solution.obj())
