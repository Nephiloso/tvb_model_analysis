from ray import tune
import numpy as np
import pandas as pd
import os
import datetime
from ray.tune.suggest.bayesopt import BayesOptSearch
from one_region_simu import *
from evaluation import *
import pdb
import sys
from CustomLoggerCallback import *

data_folder = os.path.join(os.getcwd(), 'data')
exp_name = '20220623_trial1'
ray_folder = os.path.join(data_folder,'ray_results')
exp_folder = os.path.join(ray_folder, exp_name)
restore_folder = exp_folder
os.makedirs(exp_folder, exist_ok=True)

config={
    'nsig':tune.uniform(0.0001,0.004),
    'sigma': tune.uniform(0.2, 0.5),
        "c_ee": tune.uniform(8, 18),
        'c_ie': tune.uniform(7.5, 10)}

def run_model(config):
    params = config_params({"sigma":config['sigma'], "c_ee":config['c_ee'],'c_ie':config['c_ie']})
    # score = eva_limit_cycle(params)
    # if not score:
    region = config_one_region()
    surface = config_surface(region)
    sim = config_simulator(params, region, surface)
    result_name = run_simulation(sim, params, check_point=12500,data_folder=data_folder)
    score, dfa_all, f_peak = dfa_analysis(result_name)
    tune.report(score=score, dfa_all = dfa_all, f_peak=f_peak)
    
smoke_test = False
restore = True
nsample = -1
debug = False
bayesopt = BayesOptSearch(metric="score",
                          mode="min",
                          points_to_evaluate=[{"nsig":0.0001,"sigma": 0.3,"c_ee": 8, 'c_ie':8.5}],
                          random_search_steps=150)
# pdb.set_trace()
# if restore:
#     bayesopt.restore(os.path.join(restore_folder, 'baye_checkpoint.pkl'))
if debug:
    sys.path.append('/opt/app-root/src/.local/lib/python3.8/site-packages/ray/tune')
    pdb.set_trace()
analysis = tune.run(run_model,
                    verbose=3,
                    name = exp_name,
                    local_dir=ray_folder,
                    search_alg=bayesopt,
                    # config=config,
                    metric="score",
                    mode="min",
                    progress_reporter=tune.CLIReporter(metric_columns='score',parameter_columns=config),
                    callbacks=[CustomLoggerCallback(bayesopt=bayesopt, filefolder=exp_folder, filename= "log_test.txt")],
                    stop={'score':0.06},
                    # resume= 'ERRORED_ONLY', # do this when you whant to resume a stopped experiment option:'ERRORED_ONLY'
                    resume= True,
                    # max_concurrent_trials = 3,
                    num_samples=2 if smoke_test else nsample)