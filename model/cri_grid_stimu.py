import numpy as np
import pandas as pd
import os
import datetime
from one_region_simu import *
from evaluation import *
import sys
from ray import tune
import dfa
'''
applying stimulus to the model with a variety # of nodes for the whole parameter space.
This is the first half of the param space.
'''
data_folder = '/scratch/grid_stimu'
config={"stim_size": tune.grid_search([1, 10, 40, 100, 120, 150, 170]),"c_ee": tune.grid_search([i for i in range(6,23,4)]),
         "c_ei": tune.grid_search([i for i in range(6,15,4)])}
def run_model(config):
    params = config_params(comb=['c_ee','c_ei','stim_size'],nsig=0.0004, stim_size=config['stim_size'], sigma=1,c_ee=config['c_ee'],c_ei=config['c_ei'],
                        c_ie= 9.477017434965216,c_ii=7.4836638161762013,a_e=1.3, b_e=2.8,
                        c_e=7.0,a_i=2.0,alpha_e=1)
    region = config_one_region()
    surface = config_surface(region,surface = True)
    stimulus = config_stimulus(params, surface, onset=750, mode= 'RandomPulses', amp=0.001, temp_path=data_folder)
    sim = config_simulator(params, region, surface, integ_mode= 'stochastic', stimulus=stimulus)
    result_name = run_simulation(sim, params,data_folder=data_folder)
    score, dfa_all, f_peak = dfa_analysis(result_name, not_remove=True)
    tune.report(params=params, score=score, dfa_all = dfa_all, f_peak=f_peak)

analysis = tune.run(run_model,
                        verbose=3,
                        name = 'grid_stimu',
                        local_dir=data_folder,
                        config=config,
                        metric="score",
                        progress_reporter=tune.CLIReporter(metric_columns='score',parameter_columns=config),
                        max_concurrent_trials = 3)
