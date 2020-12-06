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
                                                                '(just_const, just_rconst, grasp, lsearch, vnd, sa, gvns)')
    parser.add_argument("--inst_file", type=str, default=default_inst_file,
                        help='problem instance file')
    parser.add_argument("--incremental_delta", type=str, default="True", help="")
    parser.add_argument("--neighborhood", type=str, default="2opt", help="for local? (2opt, 3opt, 2.5opt, xchg, insert, sblock...)")
    parser.add_argument("--step", type=str, default="best", help="(best, first, random) improvement")
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
    
    feasible = []
    feasible_times = []
        
    while elapsed < settings.mh_ttime:
        iter_start = time.process_time()
        solution.construct(Construct.GREEDY_EDGE_RANDOM, alpha=alpha_val)
        iter_end = time.process_time()
        elapsed = iter_end - start
#        print(elapsed, "alpha", alpha_val, "obj", solution.obj(), "best", best_sol.obj())

        obj_val = abs(solution.obj())
        
        if solution.obj() <= solution.inst.valid_threshold:
            feasible.append(obj_val)
            feasible_times.append(iter_end - iter_start)
        
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
                
        if solution.is_better(best_sol):
#            print("success at", alpha_val)
            best_sol, solution = solution, best_sol
            
        alpha_val = random.uniform(min_alpha, min(1.0, avg_alpha + max(min_alpha, 2*avg_alpha)) )
            
        iterations += 1
        
    count = len(feasible)
    obj_std = np.std(feasible) if count else 0
    obj_mean = np.mean(feasible) if count else 0
    obj_min = min(feasible) if count else 0
    obj_max = max(feasible) if count else 0
    t_std = np.std(feasible_times) if count else 0
    t_mean = np.mean(feasible_times) if count else 0
    t_min = min(feasible_times) if count else 0
    t_max = max(feasible_times) if count else 0
    
    times = (t_mean, t_min, t_max, t_std)
    objs = (obj_mean, obj_min, obj_max, obj_std)
    
    return elapsed, iterations, count, objs, times, best_sol

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
        total_time, iterations, num_feasible, objs, times, best_sol = xrun_rconst(solution, settings.mh_ttime)
        
        s = f"rconst;na;{len(best_sol.x)};{total_time};{iterations};{num_feasible};"
        obj_mean, obj_min, obj_max, obj_std = objs
        s += f"{obj_mean};{obj_min};{obj_max};{obj_std};"
        
        t_mean, t_min, t_max, t_std = times
        s += f"{t_mean};{t_max};{t_min};{t_std};"
        
        s += f"{abs(best_sol.obj())};"
        s += " ".join(str(v) for v in best_sol.x)
        
        print("~~~solution~~~", s)       
        
        best_sol.check()
    elif settings.alg == 'grasp':
        alg = GRASP(solution, Method("rconst", CBTSPSolution.random_construct, {'alpha': 0.0} ), 
                    Method(f"li_{settings.neighborhood}_{settings.step}", CBTSPSolution.local_improve, 
                           (neighborhood.incremental_delta(incremental_eval), step_func)))
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
        
        print("~~~solution~~~", format_solution("grasp", f"{settings.neighborhood}_{settings.step}", alg.run_time, alg.iteration, alg.incumbent))
        
    elif settings.alg.startswith('lsearch'):
        #if no shake method is given it always stops after one non-improving iteration
        #no matter what. So using a dummy to make it work for random step, but not a great solution.
        #There is also no way to keep searching past a local optimum with other step
        #functions as no worse solution than current best is ever tried tried and they always
        #return the same solution in a local optimum.
        ls_settings = {'mh_tciter': -1} 
        
        li_meth = [Method(f"li_{settings.neighborhood}_{settings.step}", CBTSPSolution.local_improve, 
                           (neighborhood.incremental_delta(incremental_eval), step_func))]
        
        if settings.alg == 'lsearch-r':
            sh_meth = [Method(f"sh_{settings.neighborhood}_random", CBTSPSolution.local_improve, 
                           (neighborhood.incremental_delta(incremental_eval), Step.RANDOM))]
        else:
            sh_meth = []
        
        alg = GVNS(solution, [Method("ch_ham_path", CBTSPSolution.construct, Construct.HAMILTON_PATH)],
                   li_meth, sh_meth, ls_settings)

        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
        
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
#            'mh_tciter': 10000, # Shortcut: Abort after 10000 non-improving iterations - remove for real tests

#             should be set in command line probably
#            'mh_ttime': 15*60 # Limited to 15 min CPU time
        }
        alg = SA_CBTSP(solution, [Method("rconst", CBTSPSolution.construct, Construct.GREEDY_EDGE)], CBTSPSolution.random_move_delta_eval, CBTSPSolution.apply_neighborhood_move, None, sa_settings)
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()
    elif settings.alg == "vnd":
        # Which neighborhoods do we want to use here?

        # 3-opt is too slow on to use on larger instances, but could maybe only enable it on smaller ones?
        # 2.5-opt and insertion together does not make sense as it's just 2-opt plus some insertion moves
        # Otherwise maybe restrict to : 2opt, insert, xchg, short block ?
        ms = [ Method("li_2opt_best", CBTSPSolution.local_improve, (NeighborhoodSpec.TWO_OPT, Step.BEST)),
#              Method("li_2.5opt_best", CBTSPSolution.local_improve, (NeighborhoodSpec.TWO_HALF_OPT, Step.BEST)),
               Method("li_xchg_best", CBTSPSolution.local_improve, (NeighborhoodSpec.TWO_XCHG, Step.BEST)),
               Method("li_smove_best", CBTSPSolution.local_improve, (NeighborhoodSpec.SINGLE_INSERT, Step.BEST)),
               Method("li_sblock_best", CBTSPSolution.local_improve, (NeighborhoodSpec.SHORT_BLOCK, Step.BEST)) ]
        
        # at 400 points one iteration of 3opt best neighbor search takes almost a minute on my PC
        threshold = 500 
        if len(solution.x) < threshold:
            ms += [Method("li_3opt_best", CBTSPSolution.local_improve, (NeighborhoodSpec.THREE_OPT, Step.BEST))]
        
        # Not sure about construction method, if run repeatedly random construction can sometimes give better 
        # results with small instances but larger instances will be bad or infeasible
        alg = GVNS(solution, [Method("rconst", CBTSPSolution.construct, Construct.GREEDY_EDGE_RANDOM)], random.sample(ms, len(ms)), [])
        alg.run()
        logger.info("")
        alg.method_statistics()
        alg.main_results()

#    print(solution, solution.obj())
