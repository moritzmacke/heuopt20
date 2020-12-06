import logging
from typing import Callable
import time

from distutils.util import strtobool
from sortedcontainers import SortedList, SortedSet, SortedDict

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
    parser.add_argument("--incremental_delta", type=str, default="True", help="")
    parser.add_argument("--neighborhood", type=str, default="2opt", help="for local? (2opt, 3opt, 2.5opt, xchg, insert, sblock...)")
    parser.add_argument("--step", type=str, default="best", help="(best, first, random) improvement")
    parser.add_argument("--out_file", type=str, default="", help="file to write final result")
    parse_settings(seed=seed)


param_to_neighbor = {'2opt': NeighborhoodSpec.TWO_OPT,
                     '3opt': NeighborhoodSpec.THREE_OPT,
                     '2.5opt': NeighborhoodSpec.TWO_HALF_OPT,
                     'xchg': NeighborhoodSpec.TWO_XCHG,
                     'insert': NeighborhoodSpec.SINGLE_INSERT,
                     'sblock': NeighborhoodSpec.SHORT_BLOCK }

param_to_step = {'best': Step.BEST,
                 'first': Step.FIRST,
                 'random': Step.RANDOM }

def format_solution(method, params, total_time, iterations, solution):
    s = f"{method};{params};{len(solution.x)};{abs(solution.obj())};{total_time};{iterations};"
    s += " ".join(str(v) for v in  solution.x)
    return s

def run_rconst(solution, timeout):
    best_sol = solution.copy()
    start = time.process_time()
    elapsed = 0
        
    alpha_val = 0
    alpha_step = 1.0/len(solution.inst.edges)
    min_alpha = 2*alpha_step       
    max_alpha = min_alpha
    iterations = 0
        
    while elapsed < settings.mh_ttime:
        solution.construct(Construct.GREEDY_EDGE_RANDOM, alpha=alpha_val)
        elapsed = time.process_time() - start
#        print(elapsed, "alpha", alpha_val, "obj", solution.obj(), "best", best_sol.obj())
        
        if solution.is_better(best_sol):
#            print("success at", alpha_val)
            best_sol, solution = solution, best_sol
            
        max_alpha *= 1.01#+= alpha_step
        if max_alpha > min(1.0, (alpha_step* max(100,len(solution.inst.edges)/10) )):
            max_alpha = min_alpha
        alpha_val = max_alpha#random.uniform(min_alpha, max_alpha)
        iterations += 1
    
    return elapsed, iterations, best_sol

def xrun_rconst(solution, timeout):
    best_sol = solution.copy()
    start = time.process_time()
    elapsed = 0
        
    top_solutions = SortedDict()
    
    cum_alpha = 0
        
    alpha_val = 0
    alpha_step = 1.0/len(solution.inst.edges)
    min_alpha = 2*alpha_step
    avg_alpha = alpha_step
    iterations = 0
        
    while elapsed < settings.mh_ttime:
        solution.construct(Construct.GREEDY_EDGE_RANDOM, alpha=alpha_val)
        elapsed = time.process_time() - start
#        print(elapsed, "alpha", alpha_val, "obj", solution.obj(), "best", best_sol.obj())

        obj_val = abs(solution.obj())
        if obj_val not in top_solutions and alpha_val != 0:
            cum_alpha += alpha_val
            top_solutions[obj_val] = alpha_val
            
        if len(top_solutions) > 20:
            obj_val, alpha = top_solutions.popitem()
            cum_alpha -= alpha
            avg_alpha = cum_alpha / len(top_solutions)
        else:
            avg_alpha += alpha_step
            avg_alpha = min(avg_alpha, 1.0)
        
#        print(avg_alpha, top_solutions)
#        print(avg_alpha)
        
        if solution.is_better(best_sol):
#            print("success at", alpha_val)
            best_sol, solution = solution, best_sol
            
        alpha_val = random.uniform(min_alpha, min(1.0, avg_alpha + max(min_alpha, 2*avg_alpha)) )
            
        iterations += 1
    
    return elapsed, iterations, best_sol

if __name__ == '__main__':
    
    add_general_arguments_and_parse_settings('0010.txt', 0)

    settings.mh_workers = 1 # always one thread
    settings.mh_ttime = 20 if settings.mh_ttime < 0 else settings.mh_ttime # don't allow unlimted runtime if none set...
    settings.mh_titer = -1

    init_logger()
    logger = logging.getLogger("pymhlib")
    logger.info(get_settings_as_str())
       
    instance = CBTSPInstance(inst_dir + settings.inst_file)
    logger.info("instance read:\n" + str(instance))
    solution = CBTSPSolution(instance)
    
    neighborhood = param_to_neighbor[settings.neighborhood]
    step_func = param_to_step[settings.step]
    incremental_eval = strtobool(settings.incremental_delta)
   
    if settings.alg == 'just_const':
        start = time.process_time()
        solution.construct(Construct.HAMILTON_PATH, None, timeout=settings.mh_ttime)
        total_time = time.process_time() - start
        print("~~~solution~~~",format_solution("dconst", "na", total_time, -1, solution))
        solution.check()
    elif settings.alg == 'just_rconst':
        total_time, iterations, sol = xrun_rconst(solution, settings.mh_ttime)
        print("~~~solution~~~", format_solution("rconst", "na", total_time, iterations, sol))
        sol.check()
    elif settings.alg == 'grasp':
        alg = GRASP(solution, Method("rconst", CBTSPSolution.random_construct, {'alpha': 0.0} ), 
                    Method(f"li_{settings.neighborhood}_{settings.step}", CBTSPSolution.local_improve, 
                           (neighborhood.incremental_delta(incremental_eval), step_func)))
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
        
        print("~~~solution~~~", format_solution("grasp", f"{settings.neighborhood}_{settings.step}", alg.run_time, alg.iteration, alg.incumbent))
        if settings.out_file != "":
            with open(settings.out_file, "a") as out:
                out.write("\n" + format_solution("grasp", f"{settings.neighborhood}_{settings.step}", alg.run_time, alg.iteration, alg.incumbent))
                if len(alg.incumbent.get_invalid_edges()) > 0:
                    out.write(" (!)")

    elif settings.alg == 'lsearch':
        #if no shake method is given it always stops after one non-improving iteration
        #no matter what. So using a dummy to make it work for random step, but not a great solution.
        #There is also no way to keep searching past a local optimum with other step
        #functions as no worse solution than current best is ever tried tried and they always
        #return the same solution in a local optimum.
        ls_settings = {'mh_tciter': -1 if step_func == Step.RANDOM else 1} 
        
        alg = GVNS(solution, [Method("ch_ham_path", CBTSPSolution.construct, Construct.HAMILTON_PATH)],
                   [Method(f"li_{settings.neighborhood}_{settings.step}", CBTSPSolution.local_improve, 
                           (neighborhood.incremental_delta(incremental_eval), step_func))],
                   [Method(f"sh_{settings.neighborhood}_random", CBTSPSolution.local_improve, 
                           (neighborhood.incremental_delta(incremental_eval), Step.RANDOM))], ls_settings)
        alg.run()
        logger.info("")
        alg.method_statistics()
#        alg.main_results()
        
        print("~~~solution~~~", format_solution("lsearch", f"ham_{settings.neighborhood}_{settings.step}", alg.run_time, alg.iteration, alg.incumbent))
        
    elif settings.alg == 'gvns':
        alg = GVNS(solution, [Method(f"ch0", CBTSPSolution.construct, Construct.HAMILTON_PATH)], 
                   [Method(f"li_{settings.neighborhood}_{settings.step}", CBTSPSolution.local_improve, 
                           (neighborhood.incremental_delta(incremental_eval), step_func))], 
                   [Method(f"sh{i}", CBTSPSolution.shaking, i) for i in range(1, 2)])
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
    elif settings.alg == "sa":
        sa_settings = {
            'mh_titer': -1,
            'mh_sa_T_init': solution.inst.bigM/10,
            'mh_sa_T_reheat': solution.inst.bigM/50,
            'mh_sa_reheat_iter': 5*1000*1000
        }
        alg = SA_CBTSP(solution, [Method("rconst", CBTSPSolution.construct, Construct.HAMILTON_PATH)],
                       neighborhood.incremental_delta(incremental_eval).random_move_delta, neighborhood.apply_move,
                       None, sa_settings)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()

        if settings.out_file != "":
            with open(settings.out_file, "a") as out:
                out.write("\n" + format_solution("sa",settings.neighborhood,alg.incumbent_time,alg.incumbent_iteration,alg.incumbent))
                if len(alg.incumbent.get_invalid_edges()) > 0:
                    out.write(" (!)")
    elif settings.alg == "vnd":
        # Which neighborhoods do we want to use here?

        ms = [ Method("li_2opt_best", CBTSPSolution.local_improve, (NeighborhoodSpec.TWO_OPT.incremental_delta(incremental_eval), Step.BEST)),
               Method("li_xchg_best", CBTSPSolution.local_improve, (NeighborhoodSpec.TWO_XCHG.incremental_delta(incremental_eval), Step.BEST)),
               Method("li_smove_best", CBTSPSolution.local_improve, (NeighborhoodSpec.SINGLE_INSERT.incremental_delta(incremental_eval), Step.BEST)),
               Method("li_sblock_best", CBTSPSolution.local_improve, (NeighborhoodSpec.SHORT_BLOCK.incremental_delta(incremental_eval), Step.BEST)) ]

        alg = GVNS(solution, [Method("rconst", CBTSPSolution.construct, Construct.HAMILTON_PATH)], ms, [])
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()

        if settings.out_file != "":
            with open(settings.out_file, "a") as out:
                out.write("\n" + format_solution("vnd","2opt->xchg->smove->sblock", alg.run_time, alg.iteration, alg.incumbent))
                if len(alg.incumbent.get_invalid_edges()) > 0:
                    out.write(" (!)")

#    print(solution, solution.obj())
