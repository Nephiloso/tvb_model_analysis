from one_region_simu import *
import matplotlib.pyplot as plt
import dfa
from scipy.signal import welch, hamming, find_peaks
from utils import *

BANDS = [[0.5,4],[4,8],[8,16],[16,30],[30,100]]
def eva_limit_cycle(params, dt=2**-1):
    '''Check whether the trace on phase plane forms a limit cycle with size < ~15Hz. Within 210ms, there should be more than 3 peaks in the trace, i.e. 3 periods within 210ms -> 15Hz
       Not successful. Adandoned.'''
    region = config_one_region()
    surface = config_surface(region)
    sim = config_simulator(params, region, surface, integ_mode='deterministic', simu_length=210)
    savg_time, savg_data = run_simulation(sim, params, return_signal=True)
    savg_t = np.squeeze(np.array(savg_time))
    savg = np.squeeze(np.array(savg_data))
    
    # Check number of peaks to determin whether it's periodic or monotonic
    waive_length = int(5/dt) # Leaving out the first 5ms to get rid of artifact from the random seed
    x_peak = find_peaks(savg[waive_length:])
    if len(x_peak[0])<3:
        return 2
    plot_result(savg_t,savg, peaks=x_peak[0]+waive_length)
    return None

def dfa_analysis(result_name, dt=1, not_remove = False, save_inter=False):
    df = pd.read_csv(result_name, header = None)
    data = df[1].to_numpy()
    if sum(np.isnan(data))>0:  # get rid of unrealistic trials with NaN
        os.remove(result_name)
        return 2
    R = []
    if save_inter:
        peak = get_psd_peak(data, dt=dt,save_path=result_name)
    else:
        peak = get_psd_peak(data, dt=dt)
    for i in range(len(BANDS)):
        raw = dfa.load_data([data], sfreq = 1000/dt)
        R0 , _ = dfa.compute_DFA(raw, l_freq=BANDS[i][0], h_freq=BANDS[i][1])
        R.append(R0)
        if (peak > BANDS[i][0]) & (peak < BANDS[i][1]):
            band = i
        if save_inter:
            filtered_d, filtered_t = raw[:]
            df = pd.DataFrame({'time':filtered_t, 'data':filtered_d[0]})
            filtered_name=result_name[:result_name.find('results.csv')]+'envelope'+str(BANDS[i][0])+'_'+str(BANDS[i][1])+'.csv'
            df.to_csv(filtered_name, mode='a', index=False, header=False)
        del raw
        del filtered_t
        del filtered_d
        del filtered_name
    raw = dfa.load_data([data], sfreq = 1000/dt)
    R0 , _ = dfa.compute_DFA(raw,filter_data=False)
    R.append(R0)
    dfa_all = {'delta':R[0][0],'theta':R[1][0],'alpha':R[2][0],'beta':R[3][0], 'gamma':R[4][0],'raw':R[5][0]}
    penalty0 = 0.87  # penalty for delta, gamma
    penalty1 = 0.92  # penalty for theta, beta

    if R[band] > R[5]:
        dfa0 = R[band]
    else:
        dfa0 = R[5]
    dfa_all['dfa_this_trial'] = dfa0[0]
    if (band!=2) & (band!=1) & (band!=3):
        score = abs(dfa0*penalty0-0.85)[0]
    elif (band==1) or (band==3):
        score = abs(dfa0*penalty1-0.85)[0]
    elif band==2:
        score = abs(0.85 - dfa0[0])

    if (dfa0<0.75) and not not_remove:
        os.remove(result_name)
    return (score, dfa_all, peak)

def get_score(dfa_all, peak):
    for i in range(len(BANDS)):
        if (peak > BANDS[i][0]) & (peak < BANDS[i][1]):
            band = i
    penalty0 = 0.87  # penalty for delta, beta, gamma
    penalty1 = 0.92  # penalty for theta
    R = list(dfa_all.values())
    if R[band] > R[5]:
        dfa0 = R[band]
    else:
        dfa0 = R[5]
    dfa_all['dfa_this_trial'] = dfa0
    if (band!=2) & (band!=1) & (band!=3):
        score = abs(dfa0*penalty0-0.85)
    elif (band==1) or (band==3):
        score = abs(dfa0*penalty1-0.85)
    elif band==2:
        score = abs(0.85 - dfa0[0])
    return score
    
def get_psd(data, fs):
    nfft= 2048
    overlap = 64
    f, Pxxf = welch(data, fs, window=hamming(nfft), noverlap=overlap, nfft=nfft, return_onesided=True, scaling='density')
    return (f, Pxxf)

def get_psd_peak(data, **kwargs):
    '''Get the frequency where the PSD has a peak'''
    if list(kwargs.keys())[0] == ('fs' or 'Fs'):
        fs = list(kwargs.values())[0]
    elif list(kwargs.keys())[0] == ('dt'):
        fs = 1000/list(kwargs.values())[0]
    else:
        raise Exception("Should specify frequency by providing arguement 'fs' or 'dt'.")
    f, Pxxf = get_psd(data, fs)
    idx1 = (np.abs(f - 100)).argmin()
    idx0 = (np.abs(f - 0.5)).argmin()
    if f[idx0]<0.5:
        idx0+=1
    peak_idx = np.argmax(Pxxf[idx0:idx1])
    if list(kwargs.keys())[1] == ('save_path'):
        df = pd.DataFrame({'f':f, 'Pxxf':Pxxf})
        result_name = list(kwargs.values())[1]
        psd_name=result_name[:result_name.find('results.csv')]+'psd'+'.csv'
        df.to_csv(psd_name, mode='a', index=False, header=False)
    return f[peak_idx+idx0]
