import logging
from typing import Callable

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
    parser.add_argument("--alg", type=str, default='sa', help='optimization algorithm to be used '
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
    elif settings.alg == 'just_rconst':
        solution.construct(Construct.GREEDY_EDGE_RANDOM, None)
        solution.check()
    elif settings.alg == 'grasp':
        alg = GRASP(solution, Method("rconst", CBTSPSolution.construct, Construct.GREEDY_EDGE_RANDOM), Method("search", CBTSPSolution.local_improve, (Neighbor.KOPT2, Step.BEST)), ownsettings)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
    elif settings.alg == 'lsearch':
        raise NotImplementedError
    elif settings.alg == 'gvns':
        alg = GVNS(solution, [Method(f"ch0", CBTSPSolution.construct, Construct.GREEDY_EDGE)], [Method("li_2opt_best", CBTSPSolution.local_improve, (Neighbor.KOPT2, Step.BEST))], [Method(f"sh{i}", CBTSPSolution.shaking, i) for i in range(1, 2)], ownsettings)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
    elif settings.alg == "sa":
        alg = SA_CBTSP(solution, [Method("rconst", CBTSPSolution.construct, Construct.GREEDY_EDGE_RANDOM)], CBTSPSolution.random_move_delta_eval, CBTSPSolution.apply_neighborhood_move, None)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()

#    print(solution, solution.obj())
