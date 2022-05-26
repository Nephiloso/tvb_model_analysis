from ray import tune
import numpy as np
# import pandas as pd
import os
# from datetime import date
from ray.tune.suggest.bayesopt import BayesOptSearch
from one_region_simu import *

data_folder = os.path.join(os.getcwd(), 'data')

config={ "c_ee": tune.uniform(12.5, 13.5),
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
nsample = 50 # 1000
bayesopt = BayesOptSearch(metric="score", mode="min")
analysis = tune.run(run_model,
                    name = '20220522_trail1',
                    local_dir=os.path.join(data_folder,'ray_results'),
                    search_alg=bayesopt,
                    config=config,
                    # progress_reporter=tune.JupyterNotebookReporter(score),
                    num_samples=1 if smoke_test else nsample)

print("Best hyperparameters found were: ", analysis.best_config)
# params = config_params()
# region = config_one_region()
# surface = config_surface(region)
# sim = config_simulator(params, region, surface)
# result_name, oscillation = run_simulation(sim)
# score = evaluation(result_name, oscillation)