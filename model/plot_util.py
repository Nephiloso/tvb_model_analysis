'''
Frequently used plotting functions.
'''
from scipy.signal import welch, hamming
import matplotlib.pyplot as plt
import numpy as np

def plot_psd(data, **kwargs):
    if 'nfft' in kwargs.keys():
        nfft = kwargs['nfft']
    else:
        nfft= 2048
    if nfft > data.shape[0]:
        nfft = int(data.shape[0])
    overlap = int(nfft/8)
    if ('fs' or 'Fs') in kwargs.keys():
        fs = list(kwargs.values())[0]
    elif 'dt' in kwargs.keys():
        fs = 1000/list(kwargs.values())[0]
    else:
        raise Exception("Should specify frequency by providing arguement 'fs' or 'dt'.")
    f, Pxxf = welch(data, fs, window=hamming(nfft), noverlap=overlap, nfft=nfft, return_onesided=True, scaling='density')

    plt.semilogy(f, Pxxf, '-o')
    plt.grid()
    if f[-1] >100:
        plt.xlim([1, 100])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('PSD (dB/Hz)')
    plt.show()
    return (f, Pxxf)

def plot_result(time, signal, interval=[], t_unit='ms', peaks=[]):
    #Plot region averaged time series
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(111)
    if signal.ndim == 4:
        if interval:
            ax.plot(time[interval[0]:interval[-1]], signal[interval[0]:interval[-1], 0, :, 0])
        else:
            ax.plot(time,signal[:, 0, :, 0])
    elif signal.ndim == 1:
        if interval:
            ax.plot(time[interval[0]:interval[-1]], signal[interval[0]:interval[-1]])
        else:
            ax.plot(time,signal)
    else:
        raise Exception("Dimension should be 4 or 1!")
        return None
    if np.any(peaks):
        if interval:
            x = peaks[peaks>interval[0]]
            y = peaks[peaks<interval[1]]
            np.intersect1d(x,y)
        plt.scatter(time[peaks], signal[peaks],color='r',linewidths=0.1)
    plt.title("Region average")
    plt.xlabel("time/"+t_unit)
    #Show them
    plt.show()