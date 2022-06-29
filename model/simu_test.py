from one_region_simu import *
from evaluation import *

data_folder = '/scratch/pwutmp0'
import os;os.makedirs(data_folder,exist_ok=True)
params = config_params(nsig=0.0010545699, sigma=0.25502, c_ee=14.99466,c_ie=7.185087)
region = config_one_region()
surface = config_surface(region)
sim = config_simulator(params, region, surface)
result_name = run_simulation(sim, params,data_folder=data_folder)
score, dfa_all, f_peak = dfa_analysis(result_name,not_remove = True)
print(score)
print(dfa_all)
print(f_peak)
