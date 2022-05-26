from tvb.simulator.lab import *
import pandas as pd
from datetime import date
import os
import dfa
import numpy as np

'''
This script is a neural mass model for V1.
'''

def config_params(params=[]):
    # default settings
    k_e = np.array([1.0])            #Max value of excitatory response function
    k_i = np.array([2.0])            #Max value of inhibitory response function
    r_e = np.array([1.0])            #Excitatory refractory period
    r_i = np.array([1.0])            #Inhibitory refractory period
    c_ee = np.array([13.25])          #Excitatory to excitatory coupling coefficient
    c_ei = np.array([12.0])          #Inhibitory to excitatory coupling coefficient
    c_ie = np.array([9.23])         #Excitatory to inhibitory coupling coefficient
    c_ii = np.array([2.0])          #Inhibitory to inhibitory coupling coefficient
    tau_e = np.array([10.0])         #Membrane time-constant for the excitatory population
    tau_i = np.array([10.0])         #Membrane time-constant for the inhibitory population
    a_e = np.array([1.3])            #Slope parameter for the excitatory response function
    b_e = np.array([2.8])           #Position of the maximum slope of the excitatory sigmoid function
    c_e = np.array([7.0])            #Amplitude parameter for the excitatory response function
    a_i = np.array([2.0])            #Slope parameter for the inhibitory response function
    b_i = np.array([4.0])            #Position of the maximum slope of the inhibitory sigmoid function
    c_i = np.array([1.0])            #Amplitude parameter for the inhibitory response function
    theta_e = np.array([2.0])        #Excitatory treshold
    theta_i = np.array([1.5])        #Inhibitory treshold
    alpha_e = np.array([1.0])        # ?Balance parameter between excitatory and inhibitory masses
    alpha_i = np.array([1.0])        # ?Balance parameter between excitatory and inhibitory masses
    P = np.array([2.5])              #External stimulus to the excitatory population
    Q = np.array([0.0])              #External stimulus to the inhibitory population
    shift_sigmoid = np.array([False])
                         
    DEFAULT_PARAMS = {'k_e':k_e, 'k_i':k_i, 'r_e':r_e, 'r_i':r_i, 'c_ee':c_ee, 'c_ei':c_ei,
                      'c_ie':c_ie, 'c_ii':c_ii, 'tau_e':tau_e, 'tau_i':tau_i,
                      'a_e':a_e, 'b_e':b_e, 'c_e':c_e, 'a_i':a_i, 'b_i':b_i, 'c_i':c_i,
                      'theta_e':theta_e, 'theta_i':theta_i, 'alpha_e':alpha_e, 'alpha_i':alpha_i,
                      'P':P, 'Q':Q, 'shift_sigmoid':shift_sigmoid}
    # update the params
    if params:
        for i,j in enumerate(params.keys()):
            if j in DEFAULT_PARAMS:
                DEFAULT_PARAMS[j] = np.array([params[j]])

    return DEFAULT_PARAMS

def config_one_region():
    # region-based modeling
    white_matter = connectivity.Connectivity.from_file()
    tract_lengths = white_matter.tract_lengths[73:74,73:74]
    region_labels = white_matter.region_labels[73:74]
    centres = white_matter.centres[73:74]
    cortical = white_matter.cortical[73:74]
    orientations = white_matter.orientations[73:74]
    areas = white_matter.areas[73:74]
    weights = white_matter.weights[73:74,73:74]
    one_region = connectivity.Connectivity(number_of_regions = 1,
                                            tract_lengths = tract_lengths,
                                            region_labels = region_labels,
                                            centres = centres,
                                            cortical = cortical,
                                            orientations = orientations,
                                            areas = areas,weights = weights)
    one_region.configure()
    return one_region


def config_surface(region, folder_path='/mnt/user/drive/My Libraries/tutorials&explorations/data/cortex/', source_file='cortex_180.zip', mapping_file='regionMapping_180_1.txt', local_file='local_connectivity_180.mat'):
    source_path = os.path.join(folder_path, source_file)
    mapping_path = os.path.join(folder_path, mapping_file)
    local_connectivity_path = os.path.join(folder_path, local_file)

    #Initialise a surface
    default_cortex = cortex.Cortex.from_file(source_file=source_path, region_mapping_file=mapping_path)
    default_cortex.coupling_strength = np.array([2**-10])
    default_cortex.local_connectivity = local_connectivity.LocalConnectivity.from_file(local_connectivity_path)
    default_cortex.region_mapping_data.connectivity = region
    default_cortex.configure()
    return default_cortex

def config_simulator(params, region, surface):
    simulator.Simulator()
    # specify the coupling
    c=coupling.Sigmoidal()
    # The algorithm to solve differential equations
    integ = integrators.HeunStochastic(dt=0.5) # a. simulate the noise of eeg recording; b. Heun is a more accurate approximation than Euler

    # specify the model
    mod = models.WilsonCowan(k_e = params['k_e'],            #Max value of excitatory response function
                         k_i = params['k_i'],            #Max value of inhibitory response function
                         r_e = params['r_e'],            #Excitatory refractory period
                         r_i = params['r_i'],            #Inhibitory refractory period
                         c_ee = params['c_ee'],          #Excitatory to excitatory coupling coefficient
                         c_ei = params['c_ei'],          #Inhibitory to excitatory coupling coefficient
                         c_ie = params['c_ie'],         #Excitatory to inhibitory coupling coefficient
                         c_ii = params['c_ii'],          #Inhibitory to inhibitory coupling coefficient
                         tau_e = params['tau_e'],         #Membrane tie-constant for the excitatory population
                         tau_i = params['tau_i'],         #Membrane time-constant for the inhibitory population
                         a_e = params['a_e'],            #Slope parameter for the excitatory response function
                         b_e = params['a_i'],           #Position of the maximum slope of the excitatory sigmoid function
                         c_e = params['c_e'],            #Amplitude parameter for the excitatory response function
                         a_i = params['a_i'],            #Slope parameter for the inhibitory response function
                         b_i = params['b_i'],            #Position of the maximum slope of the inhibitory sigmoid function
                         c_i = params['c_i'],            #Amplitude parameter for the inhibitory response function
                         theta_e = params['theta_e'],        #Excitatory treshold
                         theta_i = params['theta_i'],        #Inhibitory treshold
                         alpha_e = params['alpha_e'],        # ?Balance parameter between excitatory and inhibitory masses
                         alpha_i = params['alpha_i'],        # ?Balance parameter between excitatory and inhibitory masses
                         P = params['P'],              #External stimulus to the excitatory population
                         Q = params['Q'],              #External stimulus to the inhibitory population
                         shift_sigmoid = params['shift_sigmoid']
                        )
    #Initialise some Monitors with period in physical time
    # TODO configure monitor for each node after optimization is done
    mon_savg = monitors.SpatialAverage(period=2**-1)
    sim = simulator.Simulator(model = mod, connectivity = region,
                          coupling = c, 
                          integrator = integ, monitors = (mon_savg,),
                          simulation_length = 1000000,
                          surface = surface)
    sim.configure()
    sim.integrator.noise.nsig = np.array([1.]) # configure the noise otherwise ERROR: Bad Simulator.integrator.noise.nsig shape
    return sim

def gen_simu_name(data_folder, params):
    '''Name the file systemetically'''
    # TODO make a decorator?
    today = date.today().isoformat()
    simu_name = today +'_c_ee'+str(params['c_ee'][0])+'_c_ei'+str(params['c_ei'][0])
    result_name = os.path.join(data_folder,simu_name+"_results.csv")
    if not check_simu_name(data_folder, result_name):
        return result_name
    else:
        simu_name = today +'_c_ee'+str(params['c_ee'][0])+'_c_ei'+str(params['c_ei'][0])+'_c_ie'+str(params['c_ie'][0])
        result_name = os.path.join(data_folder,simu_name+"_results.csv")
        if not check_simu_name(data_folder, result_name):
            return result_name
        else:
            simu_name = today +'_c_ee'+str(params['c_ee'][0])+'_c_ei'+str(params['c_ei'][0])+'_c_ie'+str(params['c_ie'][0])+'_c_ii'+str(params['c_ie'][0])
            result_name = os.path.join(data_folder,simu_name+"_results.csv")
            if not check_simu_name(data_folder, result_name):
                return result_name
            else:
                simu_name = today +'_c_ee'+str(params['c_ee'][0])+'_c_ei'+str(params['c_ei'][0])+'_c_ie'+str(params['c_ie'][0])+'_c_ii'+str(params['c_ie'][0])+'_tau_e'+str(params['tau_e'][0])
                result_name = os.path.join(data_folder,simu_name+"_results.csv")
                if not check_simu_name(data_folder, result_name):
                    return result_name
                else:
                    from datetime import datetime
                    now = datetime.now()
                    current_time = now.strftime("%H:%M")
                    simu_name = today +'_c_ee'+str(params['c_ee'][0])+'_c_ei'+str(params['c_ei'][0])+'_c_ie'+str(params['c_ie'][0])+'_c_ii'+str(params['c_ie'][0])+'_tau_e'+str(params['tau_e'][0]+current_time)
                    result_name = os.path.join(data_folder,simu_name+"_results.csv")
                    return result_name
    
def check_simu_name(data_folder, result_name):
    for _, _, files in os.walk(data_folder):
        for file in files:
            if file.__contains__(result_name):
                return True
            else:
                return False

def run_simulation(sim, data_folder, params):
    savg_data = []
    savg_time = []
    result_name = gen_simu_name(data_folder, params)

    batch = 0
    check_point = 20000
    oscillation = False
    for savg in sim():
        if not savg is None:
            savg_time.append(savg[0][0])
            savg_data.append(savg[0][1])
            if len(savg_time) == check_point:
                SAVG = np.array(savg_data)
                v1 = SAVG[:,0,0,0]
                # v2 = SAVG[:,0,1,0]
                df = pd.DataFrame({'time':savg_time,
                                   's_v1':v1})
                if batch == 0:
                    df.to_csv(result_name, index=False, header=False)
                elif batch == 1:
                    if sum(abs(v1[-10000:])) < 200:
                        print("Oscillation is not generated. Trial will be terminated.")
                        sys.stdout.flush()
                        break
                    else:
                        oscillation = True
                        df.to_csv(result_name, mode='a', index=False, header=False)
                else:
                    df.to_csv(result_name, mode='a', index=False, header=False)
                batch +=1
                df= []
                savg_data = []
                savg_time = []
    return (result_name, oscillation)

def evaluation(result_name, oscillation):
    if oscillation:
        df = pd.read_csv(result_name, header = None)
        data = df[1].to_numpy()
        raw = dfa.load_data([data])
        R , _ = dfa.compute_DFA(raw)
        if R>0.7 and R<1:
            print('Criticality! R=%d'%R, flush=True)
            sys.stdout.flush()
        else:
            os.remove(result_name)
        return abs(R-0.85)
    else:
        return 10