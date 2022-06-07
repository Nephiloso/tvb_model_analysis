from one_region_simu import *

def eva_limit_cycle(params):
    '''Check whether the trace on phase plane forms a limit cycle by running a Deterministic Heun for xxx seconds'''
    region = config_one_region()
    surface = config_surface(region)
    sim = config_simulator(params, region, surface, integ_mode='deterministic', simu_length=2500)
    savg_time, savg_data = run_simulation(sim, params, return_signal=True)
    savg_t = np.squeeze(np.array(savg_time))
    savg = np.squeeze(np.array(savg_data))
    # Check numerical derivitive to the curve
    num_diri = np.diff(savg[-50:-1])/0.5
    if np.mean(num_diri)<0.1:
        return 10
    else:
        return None
    # result_name, oscillation = run_simulation(sim, data_folder, params,check_point = 12500)
    # score = evaluation(result_name, oscillation)
    
def score(result_name):
    df = pd.read_csv(result_name, header = None)
    data = df[1].to_numpy()
    raw = dfa.load_data([data])
    R , _ = dfa.compute_DFA(raw)
    if R>0.7 and R<1:
        print('Criticality! R=%d'%R, flush=True)
        sys.stdout.flush()
    # else:
    #     os.remove(result_name)
    return abs(R-0.85)[0]