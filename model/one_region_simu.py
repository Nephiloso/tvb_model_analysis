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
    '''
    # TODO Change the function as same form as follows:
    def default_pars(**kwargs):
      pars = {}

      # Excitatory parameters
      pars['tau_E'] = 1.     # Timescale of the E population [ms]
      pars['a_E'] = 1.2      # Gain of the E population
    ...
      pars['rI_init'] = 0.2  # Initial value of I

      # External parameters if any
      for k in kwargs:
          pars[k] = kwargs[k]

      # Vector of discretized time points [ms]
      pars['range_t'] = np.arange(0, pars['T'], pars['dt'])

      return pars
    '''
    # default settings
    ntau = 0 # 0.5?
    sigma = 0.3
    nsig = np.array([0.0004])
    k_e = np.array([1.0])            #Max value of excitatory response function
    k_i = np.array([2])            #Max value of inhibitory response function
    r_e = np.array([1.0])            #Excitatory refractory period
    r_i = np.array([1.0])            #Inhibitory refractory period
    c_ee = np.array([13.25])          #Excitatory to excitatory coupling coefficient
    c_ei = np.array([12.5])          #Inhibitory to excitatory coupling coefficient
    c_ie = np.array([9.23])         #Excitatory to inhibitory coupling coefficient
    c_ii = np.array([2.0])          #Inhibitory to inhibitory coupling coefficient
    tau_e = np.array([10.0])         #Membrane time-constant for the excitatory population
    tau_i = np.array([10.0])         #Membrane time-constant for the inhibitory population
    a_e = np.array([1.2])            #Slope parameter for the excitatory response function
    b_e = np.array([2.2])           #Position of the maximum slope of the excitatory sigmoid function
    c_e = np.array([4.0])            #Amplitude parameter for the excitatory response function
    a_i = np.array([0.8])            #Slope parameter for the inhibitory response function
    b_i = np.array([4.0])            #Position of the maximum slope of the inhibitory sigmoid function
    c_i = np.array([1.0])            #Amplitude parameter for the inhibitory response function
    theta_e = np.array([2.0])        #Excitatory treshold
    theta_i = np.array([1.5])        #Inhibitory treshold
    alpha_e = np.array([0.55])        # ?Balance parameter between excitatory and inhibitory masses
    alpha_i = np.array([1.0])        # ?Balance parameter between excitatory and inhibitory masses
    P = np.array([2.5])              #External stimulus to the excitatory population
    Q = np.array([0.0])              #External stimulus to the inhibitory population
    shift_sigmoid = np.array([False])

    DEFAULT_PARAMS = {'ntau': ntau, 'sigma': sigma, 'nsig': nsig, 'k_e':k_e, 'k_i':k_i, 'r_e':r_e, 'r_i':r_i, 'c_ee':c_ee, 'c_ei':c_ei,
                      'c_ie':c_ie, 'c_ii':c_ii, 'tau_e':tau_e, 'tau_i':tau_i,
                      'a_e':a_e, 'b_e':b_e, 'c_e':c_e, 'a_i':a_i, 'b_i':b_i, 'c_i':c_i,
                      'theta_e':theta_e, 'theta_i':theta_i, 'alpha_e':alpha_e, 'alpha_i':alpha_i,
                      'P':P, 'Q':Q, 'shift_sigmoid':shift_sigmoid}
    # update the params
    if params:
        for i,j in enumerate(params.keys()):
            if j in DEFAULT_PARAMS:
                if j == 'ntau':
                    DEFAULT_PARAMS[j] = params[j]
                elif j == 'sigma':
                    DEFAULT_PARAMS[j] = params[j]
                else:
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


def config_surface(region, surface = False, folder_path='/mnt/user/drive/My Libraries/tutorials&explorations/data/cortex/', source_file='cortex_180.zip', mapping_file='regionMapping_180_1.txt', local_file='local_connectivity_180.mat'):
    '''
    Configure the V1 cortex.
    # TODO Change the name of function into config_cortex and and another function to configure the surface!
    '''
    source_path = os.path.join(folder_path, source_file)
    mapping_path = os.path.join(folder_path, mapping_file)
    local_connectivity_path = os.path.join(folder_path, local_file)
    
    default_cortex = cortex.Cortex.from_file(source_file=source_path, region_mapping_file=mapping_path)
    default_cortex.coupling_strength = np.array([2**-10])
    default_cortex.local_connectivity = local_connectivity.LocalConnectivity.from_file(local_connectivity_path)
    default_cortex.region_mapping_data.connectivity = region
    if surface:
            #Initialise a surface
        surf = surfaces.CorticalSurface.from_file(source_file=source_path)
        surf.configure()
        default_cortex.region_mapping_data.surface = surf
    default_cortex.configure()
    return default_cortex

def config_model(params):
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
    return mod

def config_stimulus(cortex, stim_size, onset=10, alpha=0.5, beta=2):
    import random
    eqn_t = equations.Alpha()  # belongs to a family of exponential function, used for visual stimulus # TODO Literature search!
    eqn_t.parameters['onset'] = onset  # ! Time point of stimulus onset!
    eqn_t.parameters['alpha'] = alpha
    eqn_t.parameters['beta'] = beta
    
    eqn_s = equations.DiscreteEquation()
    
    random.shuffle(deck)
    focal = np.array(deck[:stim_size])
    stimulus = patterns.StimuliSurface(surface=cortex.surface,
                                       focal_points_triangles=focal,
                                       temporal=eqn_t,
                                       spatial = eqn_s)
    return stimulus

def config_simulator(params, region, surface, integ_mode='stochastic', simu_length=1000000, dt=1, stimulus=None):
    simulator.Simulator()
    # specify the coupling
    c=coupling.Sigmoidal()
    
    if integ_mode == 'stochastic':
        # nois = noise.Multiplicative(ntau=params['ntau'], nsig = params['nsig'], b=equations.Gaussian())  # ? Do we need to optimize ntau?
        # nois.configure_coloured(dt=0.1, shape=1)  # ? shape?? 
        # nois = noise.Multiplicative(nsig = params['nsig'], b=equations.Gaussian())
        nois = noise.Multiplicative(nsig = params['nsig'], b=equations.Gaussian(parameters={"amp": 1.0,
                                                                                      "sigma": params['sigma'],
                                                                                      "midpoint": 0.0,
                                                                                      "offset": 1.0}))
        integ = integrators.HeunStochastic(dt=0.1, noise=nois)
    elif integ_mode == 'deterministic':
        integ = integrators.HeunDeterministic(dt=0.1)
    else:
        raise ValueError("`integ_mode` should be either 'stochastic' or 'deterministic'!")

    mod = config_model(params)
    #Initialise some Monitors with period in physical time
    # TODO configure monitor for each node after optimization is done
    mon_savg = monitors.SpatialAverage(period=dt)
    if stimulus is not None:
        sim = simulator.Simulator(model = mod, connectivity = region, coupling = c,integrator = integ,
                                  monitors = (mon_savg,), simulation_length = simu_length,
                                  stimulus = stimulus, surface = surface)
    else:
        sim = simulator.Simulator(model = mod, connectivity = region, coupling = c, integrator = integ,
                                  monitors = (mon_savg,), simulation_length = simu_length,
                                  surface = surface)
    sim.configure()
    # sim.integrator.noise.nsig = params['nsig'] # configure the noise seperately otherwise ERROR: Bad Simulator.integrator.noise.nsig shape
    # if sim.log ....Bad Simulator.integrator.noise.nsig shape...
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
    return False

def run_simulation(sim, params, check_point=20000, return_signal=False, data_folder=''):
    savg_data = []
    savg_time = []

    batch = 0
    if not return_signal:
        result_name = gen_simu_name(data_folder, params)
        if not check_file_exist(result_name):
            for savg in sim():
                if not savg is None:
                    savg_time.append(savg[0][0])
                    savg_data.append(savg[0][1])
                    if len(savg_time) == check_point:
                        SAVG = np.array(savg_data)
                        v1 = SAVG[:,0,0,0]
                        # v2 = SAVG[:,0,1,0]
                        if np.isnan(v1[-1]):
                            return result_name
                        df = pd.DataFrame({'time':savg_time,
                                           's_v1':v1})
                        if batch == 0:
                            df.to_csv(result_name, index=False, header=False)
                        else:
                            df.to_csv(result_name, mode='a', index=False, header=False)
                        batch +=1
                        df= []
                        savg_data = []
                        savg_time = []
        return result_name
    else:
        for savg in sim():
            if not savg is None:
                savg_time.append(savg[0][0])
                savg_data.append(savg[0][1])
        savg_data = np.squeeze(savg_data)
        savg_time = np.squeeze(savg_time)
        return (savg_time, savg_data)
    
def check_file_exist(result_name, dt = 1, simu_length=1000000):
    if os.path.exists(result_name):
        df = pd.read_csv(result_name, header = None)
        if len(df) == simu_length/dt:
            return True
        else:
            os.remove(result_name)
            return False