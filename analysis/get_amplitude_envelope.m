function [AmplitudeEnvelope] = get_amplitude_envelope(data, fs, hp, lp)
    filter_order = 2/hp;
    AmplitudeEnvelope = abs(hilbert(filter_fir(data(:,:),hp,lp,fs,filter_order)));
end