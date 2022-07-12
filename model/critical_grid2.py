import numpy as np
import pandas as pd
import os
import datetime
from one_region_simu import *
from evaluation import *
import sys
from ray import tune
import dfa

data_folder = os.path.join(os.getcwd(),'/scratch/critical_states2')
config={
        "c_ee": tune.grid_search([i for i in range(6,23,4)]),
        "c_ei": tune.grid_search([i for i in range(6,23,4)])}
def run_model(config):
    params = config_params(c_ee=config['c_ee'],c_ei=config['c_ei'],nsig=0.00105457, sigma=0.255021,
                        c_ie= 7.818508667)
    region = config_one_region()
    surface = config_surface(region)
    sim = config_simulator(params, region, surface)
    result_name = run_simulation(sim, params,data_folder=data_folder)
    score, dfa_all, f_peak = dfa_analysis(result_name, not_remove=True)
    tune.report(params=params, score=score, dfa_all = dfa_all, f_peak=f_peak)

analysis = tune.run(run_model,
                        verbose=3,
                        name = 'critical_states2',
                        local_dir=data_folder,
                        config=config,
                        metric="score",
                        progress_reporter=tune.CLIReporter(metric_columns='score',parameter_columns=config),
                        max_concurrent_trials = 3)
