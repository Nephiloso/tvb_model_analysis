from tvb.datatypes.equations import TemporalApplicableEquation
from tvb.basic.neotraits.api import Attr, Final
from utils import get_data_folder_path
import os
import pandas as pd
import random
import numpy as np

class RandomPulses(TemporalApplicableEquation):
    """
    Create a paulse randomly each for each time interval, offset with respect to the time axis.

    **Parameters**:
    
    * dT     :  pulse repetition period
    * onset         :  time of the first pulse
    * amp           :  amplitude of the pulse
    unit: milliseconds
    """

    equation = Final(
        label="Pulse Train",
        default="where((var>=onset)&((var-onset) < dT), amp, 0)",
        doc="""Create a paulse randomly each for each time interval""")

    parameters = Attr(
        field_type=dict,
        default=lambda: {"amp": 0.001, "onset": 750, "num": 998, "time_intv": 1000,
                         "temp_path":os.path.join(get_data_folder_path(),'test.csv')},
        label="Pulse Train Parameters")
    
    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    
    def gen_pulse(self,var):
        external_input_times=[]
        for i in range(self.parameters['num']-1):
            external_input_time = int(round(self.parameters['onset']+\
                                            random.uniform(0, 1)*(self.parameters['time_intv']/2)+self.parameters['time_intv']*i))  # ms
            if external_input_time < var[0,-1]-750:
                external_input_times.append(external_input_time)
            else:
                external_input_times = external_input_times[:i]
                break
        df = pd.DataFrame({'stim_idx':external_input_times})
        df.to_csv(self.parameters['temp_path'], index=False, header=False)
        return external_input_times
    
    def evaluate(self,var):
        '''
        Generate a discrete representation of the equation for the space
        represented by ``var``.
        '''
        pulse_idx = self.gen_pulse(var)
        _pattern = np.zeros_like(var)
        for pulse in pulse_idx:
            _pattern[0,self.find_nearest(var,pulse)] = self.parameters['amp']
        
        return _pattern