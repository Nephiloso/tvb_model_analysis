from one_region_simu import *
import matplotlib.pyplot as plt
import scipy

def eva_limit_cycle(params, Fs=2**-1):
    '''Check whether the trace on phase plane forms a limit cycle with size < ~15Hz. Within 210ms, there should be more than 3 peaks in the trace, i.e. 3 periods within 210ms -> 15Hz'''
    region = config_one_region()
    surface = config_surface(region)
    sim = config_simulator(params, region, surface, integ_mode='deterministic', simu_length=210)
    savg_time, savg_data = run_simulation(sim, params, return_signal=True)
    savg_t = np.squeeze(np.array(savg_time))
    savg = np.squeeze(np.array(savg_data))
    
    # Check number of peaks to determin whether it's periodic or monotonic
    waive_length = int(5/Fs) # Leaving out the first 5ms to get rid of artifact from the random seed
    x_peak = scipy.signal.find_peaks(savg[waive_length:])
    if len(x_peak[0])<3:
        return 2
    plot_result(savg_t,savg, peaks=x_peak[0]+waive_length)
    return None
    
def get_score(result_name):
    df = pd.read_csv(result_name, header = None)
    data = df[1].to_numpy()
    raw = dfa.load_data([data])
    if sum(np.isnan(data))>0:  # get rid of unrealistic trials with NaN
        os.remove(result_name)
        return 2
    R , _ = dfa.compute_DFA(raw)
    if R>0.7 and R<1:
        print('Criticality! R=%d'%R, flush=True)
        sys.stdout.flush()
    return abs(R-0.85)[0]

def plot_result(time, signal, interval=[], t_unit='ms', peaks=[]):
    #Plot region averaged time series
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(111)
    if interval:
        if signal.ndim == 4:
            ax.plot(time[interval[0]:interval[-1]], signal[interval[0]:interval[-1], 0, :, 0])
        elif signal.ndim == 1:
            ax.plot(time[interval[0]:interval[-1]], signal[interval[0]:interval[-1]])
        else:
            print("Dimension should be 4 or 1!")
            return None
    else:
        ax.plot(time,signal)
    if np.any(peaks):
        if interval:
            x = peaks[peaks>interval[0]]
            y = peaks[peaks<interval[1]]
            np.intersect1d(x,y)
        plt.scatter(time[peaks], signal[peaks],color='r',linewidths=0.1)
    plt.title("Region average - Possible limited cycle")
    plt.xlabel("time/"+t_unit)
    #Show them
    plt.show()