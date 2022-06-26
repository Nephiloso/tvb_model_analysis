# # Customize the logger
from typing import Dict, List
import json
from ray.tune.logger import LoggerCallback
from tensorboardX import SummaryWriter
from shutil import rmtree
import os
import numpy as np
import datetime

class CustomLoggerCallback(LoggerCallback):
    """Put the results of all trails in the same file"""

    def __init__(self, filefolder, filename: str = "log.txt", bayesopt=None):
        self._filename = filename
        self._filefolder = filefolder
        self._filepath = os.path.join(filefolder, self._filename)
        self.bayesopt = bayesopt
        self.waive_params = ["time_this_iter_s","timesteps_total","episodes_total","timestamp", "pid",
                 "hostname", "node_ip", "timesteps_since_restore", "iterations_since_restore", "warmup_time"]
        
    def log_trial_start(self, trial: "Trial"):
        self._file = open(self._filepath, "at")
        self.summarywriter = SummaryWriter(os.path.join(self._filefolder,'tensor_events'))

    def log_trial_result(self, iteration: int, trial: "Trial", result: Dict):
        self._file = open(self._filepath, "at")
        for k, v in list(result.items()):
            if k in self.waive_params:
                del result[k]
                continue
            if isinstance(v, (np.ndarray,)):
                result[k] = v.tolist()
        self._file.write(json.dumps(result))
        self.summarywriter.add_scalar('score', result['score'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('dfa_delta', result['dfa_all']['delta'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('dfa_theta', result['dfa_all']['theta'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('dfa_alpha', result['dfa_all']['alpha'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('dfa_beta', result['dfa_all']['beta'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('dfa_gamma', result['dfa_all']['gamma'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('dfa_raw', result['dfa_all']['raw'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('dfa_raw', result['dfa_all']['dfa_this_trial'], len(self.bayesopt._buffered_trial_results)+1)
        self.summarywriter.add_scalar('f_peak', result['f_peak'], len(self.bayesopt._buffered_trial_results)+1)

    def on_trial_complete(self, iteration: int, trials: List["Trial"],
                          trial: "Trial", **info):
        self._file = open(self._filepath, "at")
        self._file.write(str(datetime.datetime.now())+'\n')
        self._file.close()

        self.summarywriter.close()
