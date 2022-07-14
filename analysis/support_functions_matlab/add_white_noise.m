function [ signalWithNoise ] = add_white_noise( signalTimeSeries , sdNoise )
    noiseSignal = randn(size(signalTimeSeries,1),1)*sdNoise;
    signalWithNoise = signalTimeSeries + noiseSignal;
end

