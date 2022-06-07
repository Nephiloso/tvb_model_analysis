from ray import tune
import numpy as np
# import pandas as pd
import os
# from datetime import date
from ray.tune.suggest.bayesopt import BayesOptSearch
from one_region_simu import *
from logger import *

data_folder = os.path.join(r'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data')
exp_name = '20220531_trial1'
ray_folder = os.path.join(data_folder,'ray_results')
exp_folder = os.path.join(ray_folder, exp_name)
os.makedirs(exp_folder, exist_ok=True)

config={'ntau':tune.uniform(0.0000001,10),
        "c_ee": tune.uniform(12.5, 13.5),
        "c_ei": tune.uniform(7, 10),
        'c_ie': tune.uniform(6, 18),
        'c_ii': tune.uniform(7, 10),
        'tau_e': tune.uniform(0, 150),
        'tau_i': tune.uniform(0, 150),
        'r_e': tune.uniform(0.5, 2),
        'r_i': tune.uniform(0.5, 2)
        }

def run_model(config):
    params = config_params({"c_ee":config['c_ee'],'c_ei':config['c_ei'],'c_ie':config['c_ie'],'c_ii':config['c_ii'],
                   'tau_e':config['tau_e'],'tau_i':config['tau_i'],'r_e':config['r_e'],'r_i':config['r_i']})
    region = config_one_region()
    surface = config_surface(region)
    sim = config_simulator(params, region, surface)
    result_name, oscillation = run_simulation(sim, data_folder, params)
    score = evaluation(result_name, oscillation)
    tune.report(score)

smoke_test = True
nsample = 400
bayesopt = BayesOptSearch(metric="score", mode="min",random_search_steps=100)
restore = True
if restore:
    bayesopt.restore(os.path.join(ray_folder, 'baye_checkpoint.pkl'))
analysis = tune.run(run_model,
                    verbose=3,
                    name = exp_name,
                    local_dir=ray_folder,
                    search_alg=bayesopt,
                    # config=config,
                    metric="score",
                    mode="min",
                    callbacks=[CustomLoggerCallback(filefolder=exp_folder, filename= "log_test.txt")],
                    stop={'score':0.08},
                    # resume=True, # do this when you whant to resume a stopped experiment
                    max_concurrent_trials = 1,  # TODO Test this after the checkpoint changed
                    num_samples=2 if smoke_test else nsample)
print("Best hyperparameters found were: ", analysis.best_config)
print("Best result: ", analysis.best_result)

# save the current experiment
bayesopt.save(os.path.join(ray_folder, 'baye_checkpoint.pkl'))