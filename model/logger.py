from typing import Dict, List
import json
from ray.tune.logger import LoggerCallback
import os
import numpy as np
import datetime

class CustomLoggerCallback(LoggerCallback):
    """Put the results of all trails in the same file"""

    def __init__(self, filefolder, filename: str = "log.txt"):
        self._filename = filename
        self._filepath = os.path.join(filefolder, self._filename)
        
    def log_trial_start(self, trial: "Trial"):
        self._file = open(self._filepath, "at")
        self._file.write(json.dumps(trial.evaluated_params))

    def log_trial_result(self, iteration: int, trial: "Trial", result: Dict):
        self._file = open(self._filepath, "at")
        try:
            self._file.write(json.dumps(result))
        except TypeError:
            for idx, (k, v) in enumerate(result.items()):
                if isinstance(v, (np.ndarray,)):
                    result[k] = v.tolist()
            self._file.write(json.dumps(result))

    def on_trial_complete(self, iteration: int, trials: List["Trial"],
                          trial: "Trial", **info):
        self._file = open(self._filepath, "at")
        self._file.write(str(datetime.datetime.now())+'\n')
        self._file.close()