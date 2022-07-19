function [amplitude_bins_plf,reg_amp,phase_bins_plf,phase_reg] = get_prestim_regulation(folder_path, simu_file_name, stimu_file_name)

    [time_series,stimulus_timeseries] = get_timeseries(folder_path, simu_file_name, stimu_file_name);
%     time_series_noise = add_white_noise(time_series,3);
    time_series_noise = time_series;

    delta_lp = 4;
    delta_hp = 0.5;

    pre_stim_ms =  750;
    post_stim_ms = 750;
    poststim_point = 150;
    index_poststim_point = pre_stim_ms+poststim_point+1;%+1 to skip the stimulation time
    
    Fs=1000;

    %check that I have enough time before and after the
    %stimulus in the timeseries to cover the
    %prestim/poststim analysis
    stimulus_timeseries = stimulus_timeseries(stimulus_timeseries>=pre_stim_ms & stimulus_timeseries<=length(time_series_noise)-post_stim_ms);

    no_stims = length(stimulus_timeseries);

    %% Filter signal and extract phase and amplitude
    filtered = filter_fir(time_series_noise,delta_hp,delta_lp,Fs,2/delta_hp);

    filtered = filtered - mean(filtered);
    signal_phase = angle(hilbert(filtered));
    signal_amplitude = abs(hilbert(filtered));

    %% Split phase and amplitudes timeseries by trial
    phase_per_trial = zeros(pre_stim_ms+post_stim_ms+1, no_stims);
    amplitude_per_trial  = zeros(pre_stim_ms+post_stim_ms+1, no_stims);

    for stim_idx = 1:no_stims
        stimulation_time = floor(stimulus_timeseries(stim_idx));
        phase_per_trial(:,stim_idx)  = signal_phase(stimulation_time-pre_stim_ms:(stimulation_time + post_stim_ms))';
        amplitude_per_trial(:,stim_idx)  = signal_amplitude(stimulation_time-pre_stim_ms:(stimulation_time + post_stim_ms))';
    end

    %% Compute pre-stimulus amplitude regulation
    prestimulus_amplitude_interval = [-150,-50];

    %calculate, for all trials, the amplitude in the
    %prestimulus interval
    amplitude_prestim = zeros(size(stimulus_timeseries));
    for stim_idx=1:length(stimulus_timeseries)
        amplitude_prestim(stim_idx)=mean(signal_amplitude(stimulus_timeseries(stim_idx)+prestimulus_amplitude_interval(1):stimulus_timeseries(stim_idx)+prestimulus_amplitude_interval(2)));
    end

    %sort by prestimulus amplitude
    [~,stim_index] = sort(amplitude_prestim);

    %split into 10 percentiles, by prestimulus amplitude
    amplitude_bins_plf = zeros(10,1);                
    percentile_length = length(stim_index)/10;
    percentile_start_indexes = round(percentile_length * (0:9))+1;
    percentile_end_indexes = [percentile_start_indexes(2:10)-1 length(stim_index)];

    %compute prestimulus regulation at 150 ms poststimulus
    %compute plf for trials in each bin
    for percentile_idx=1:10
        amplitude_bins_plf(percentile_idx)=circ_r(phase_per_trial(index_poststim_point,stim_index(percentile_start_indexes(percentile_idx):percentile_end_indexes(percentile_idx)))');
    end
    reg_amp = corr((1:10)',amplitude_bins_plf,'Type','Spearman');

    %% Compute prestimulus phase dependence for all timepoints around stimulus

    %make phase bins spaced at pi/16
    phases_sort_time = -5;
    index_sort_time = pre_stim_ms+phases_sort_time;
    bins_phases = -pi:pi/16:pi;
    % each bin is assigned trials wih prestimulus phases within
    % a binwidth, centered on the bin
    bin_width = pi/4;  % ??
    for bin_idx = 1:length(bins_phases)
        orig = abs(circ_dist(phase_per_trial(index_sort_time,:),bins_phases(bin_idx)));
        pcks = orig < bin_width/2;
        noPcks(bin_idx) = nnz(pcks);
    end
    minPcks = min(noPcks);

    phase_bins_plf = zeros(length(bins_phases),1);
    for bin_idx = 1:length(bins_phases)
        orig = abs(circ_dist(phase_per_trial(index_sort_time,:),bins_phases(bin_idx)));
        pcks = orig < bin_width/2;
        pcks = find(pcks);
        pcks = pcks(1:minPcks);
        phase_bins_plf(bin_idx) = circ_r(phase_per_trial(index_poststim_point,pcks)');
    end

    phase_reg = circ_r(bins_phases',phase_bins_plf);
end