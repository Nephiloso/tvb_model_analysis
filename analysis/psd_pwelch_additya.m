nfft            = 8000 ; % for freq resolution = 0.125Hz (fs/0.125 = 8000) \ Controls freq resolution of psd (most commonly preferred resolution is 0.1Hz)
fs              = 1000;  % sampling rate of CROS data(2x10^6 timepoints = 2000s) \ Sampling rate of TVB signal

% raw data to be in Nx1 format; hamming windowed PSD (other methods also available)
[pxx,f]         = pwelch(raw_cros_data(1,:)',hamming(nfft),0,nfft,fs,'psd');

% plot
plot(f,10*log10(pxx),'k')
xlim([1 20]) % adjust for visualization
legend off
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
title('PSD-Run1-Unfiltered CROS DATA') % change according to convienience