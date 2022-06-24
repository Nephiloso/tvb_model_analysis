import mne
import nolds
import numpy as np
import random as rd
from numpy.matlib import repmat
from scipy.signal import detrend

def load_data(data, sfreq = 2000, ch_types = ['eeg'], ch_names = ['V1'], add_noise=True, noise_sd=1):
    # add noise to the data
    if add_noise:
        noise = np.random.normal(size=(len(data),))*noise_sd
        data += noise
    # convert the data into mne.Raw()
    info = mne.create_info(ch_names=ch_names, sfreq=sfreq, ch_types=ch_types)
    raw = mne.io.RawArray(data, info)
    return raw

def compute_DFA(data, filter_data = True, l_freq=8, h_freq=16,fit_interval=[5,30], compute_interval=[1,120],overlap=True,channels_to_ignore=None,method="richard",return_fitting=False):
    '''
    INPUT:
    data: An instance of mne.Raw(). Could be created by load_data().
    filter: Wether to filter the signal with FIR filter.
    methods: 'righard' or 'nolds'. 'richard' is more closed to the result of matlab script.
    '''
    if filter_data:
        # filter the data
        data.filter(l_freq, h_freq, fir_window='hamming', fir_design="firwin", verbose=0)

    # add the noise?

    # get amplitute envelop
    data.apply_hilbert(envelope=True)
    
    # calculating DFA
    try:
        signal, _ = data[:]
        num_chans, num_timepoints = np.shape(signal)
    except ValueError:
        signal = [data[:]]
        num_chans, num_timepoints = np.shape(signal)
            
    sampling_frequency = data.info['sfreq']

    if channels_to_ignore == None:
        channels_to_ignore = [False]*num_chans

    length_signal = np.shape(signal)[1]

    assert fit_interval[0]>=compute_interval[0] and fit_interval[1]<=compute_interval[1],'CalcInterval should be included in ComputeInterval'
    assert compute_interval[0]>=0.1 and  compute_interval[1]<=1000,'ComputeInterval should be between 0.1 and 1000 seconds'

    #compute DFA window sizes for the given CalcInterval
    window_sizes = np.floor(np.logspace(-1,3,40) * sampling_frequency).astype(int)#%logspace from 0.1 seccond (10^-1) to 1000 (10^3) seconds

    window_sizes = window_sizes[(window_sizes >= compute_interval[0]*sampling_frequency) & \
        (window_sizes <= compute_interval[1]*sampling_frequency)]

    #get the positions of the first and last window sizes used for fitting
    fit_interval_first_window = np.argwhere(window_sizes>=fit_interval[0]*sampling_frequency)[0][0]
    fit_interval_last_window = np.argwhere(window_sizes<=fit_interval[1]*sampling_frequency)[-1][0]

    dfa_array = np.zeros(num_chans)
    dfa_array[:] = np.NAN
    fluctuations = np.zeros((num_chans,len(window_sizes)))
    fluctuations[:] = np.NAN

    if max(window_sizes)<=num_timepoints:
        for ch_idx in range(num_chans):
            if channels_to_ignore[ch_idx]:
                pass

            signal_for_channel = signal[ch_idx, :]
        
            if method=="richard":
                for i_window_size in range(len(window_sizes)):
                    if overlap==True:
                        window_overlap=0.5
                    else:
                        window_overlap = 0

                    window_offset = np.floor(window_sizes[i_window_size]*(1-window_overlap))
                    all_window_index = _create_window_indices(length_signal, window_sizes[i_window_size], window_offset)
                    # First we convert the time series into a series of fluctuations y(i) around the mean.
                    demeaned_signal = signal_for_channel - np.mean(signal_for_channel)
                    # Then we integrate the above fluctuation time series ('y').
                    signal_profile = np.cumsum(demeaned_signal)

                    x_signal = signal_profile[all_window_index]
                    d_signal = detrend(x_signal)

                    # if ch_idx==0:
                    #     print([all_window_index[0,0],all_window_index[1,0]])
                    #     print(np.shape(d_signal))

                    # And compute the fluctuation around the mean as the standard deviation for each of the detrended windows, std(dSignal,1) means normalized by N instead of N-1
                    w_detrended_fluctuations = np.std(d_signal, 1)
                    # And calculate the mean fluctuation across all windows
                    fluctuations[ch_idx,i_window_size]=np.mean(w_detrended_fluctuations)

                x=np.log10(window_sizes[fit_interval_first_window:fit_interval_last_window])
                y=np.log10(fluctuations[ch_idx,fit_interval_first_window:fit_interval_last_window])
                if return_fitting:
                    model, residual, _, _, _ = np.polyfit(x, y, 1, full=True)
                else:
                    model = np.polyfit(x, y, 1)
                dfa_array[ch_idx]=model[0]

            elif method=="nolds":
                # nolds algorithm for computing dfa
                dfa_array[ch_idx],fluctuations_2d = nolds.dfa(signal_for_channel,nvals=window_sizes,overlap=overlap,debug_data=True)
                fluctuations[ch_idx,:] = fluctuations_2d[1]
                
    if (method=='richard') & return_fitting:
        return (dfa_array, fluctuations, x, y, model, residual)
    else:
        return (dfa_array, fluctuations)

def _create_window_indices(length_signal, length_window, window_offset):

    window_starts = np.arange(0,length_signal-length_window,window_offset)
    num_windows = len(window_starts)

    one_window_index = np.arange(0,length_window)
    all_window_index = repmat(one_window_index,num_windows,1).astype(int)

    all_window_index = all_window_index + repmat(np.transpose(window_starts[np.newaxis,:]),1,length_window).astype(int)

    return all_window_index