import pandas as pd
import numpy as np
import os

def read_result(result_name):
    df = pd.read_csv(result_name, header = None)
    time = df[0].to_numpy()
    data = df[1].to_numpy()
    dt = np.mean(np.diff(time[:50]))
    return (time, data, dt)

def get_result_names(folder_path):
    u = os.listdir(path=folder_path)
    result_names = []
    for file in u:
        if 'csv' in file:
            result_names.append(os.path.join(folder_path, file))
            
    return result_names

def get_data_folder_path(folder_name='data'):
    '''
    This scipt is only used on the ebrain.
    '''
    pkg_path = os.getcwd()[:os.getcwd().find('model')]
    con_path = os.path.join(pkg_path,folder_name)
    if os.path.isdir(con_path):
        return con_path
    else:
        raise FileNotFoundError("The data folder doesn't exist or has been moved to another place!")