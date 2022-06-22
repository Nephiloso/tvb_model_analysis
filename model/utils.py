import pandas as pd
import numpy as np
import os

def read_result(result_name):
    df = pd.read_csv(result_name, header = None)
    time = df[1].to_numpy()
    data = df[1].to_numpy()
    dt = np.mean(np.diff(time[:50]))*1000
    return (time, data, dt)

def get_result_names(folder_path):
    u = os.listdir(path=folder_path)
    result_names = []
    for file in u:
        if 'csv' in file:
            result_names.append(os.path.join(folder_path, file))
            
    return result_names