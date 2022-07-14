clear all
%% Run stimulus response analysis
my_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\optimiz_stim';
Stimulus_size = [1; 10; 40; 100; 120; 150; 170];
simulation_name = '2022-07-13_c_ee6_c_ei22_stim_size';
stim_sizes_name = {'1'; '10'; '40'; '100'; '120'; '150'; '170'};
amplitude_regulations{7,1}=[];
phase_regulations{7,1}=[];
plfss{7,1}=[];
amplitude_percentile_plfss{7,1}=[];
bins_phasess{7,1}=[];
phasebin_plfss{7,1}=[];
PLF = [];
Run_Number=[0;0;0;0;0;0;0]; % 0 for noncritical network
% Run_Number=[1;1;1;1;1;1;1]; % 1 for critical network
for i=1:length(stim_sizes_name)
    simulation_data_file_path = fullfile(my_path,[simulation_name,stim_sizes_name{i},'_stim_results.csv']);  
    %a. read a file name from the folder; compare the filename and find all the
    %files from the same parameter sets
    %b. 

    % add the analysis folder to path
    addpath(genpath('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\Arthur 2020 code\model code\data_analysis'))
    result = readmatrix(simulation_data_file_path);
    v1 = result(:,2);
    raw_signal_to_process = v1;
    delta_lp = 4;
    delta_hp = 0.5;
    % result = readmatrix('C:\Users\wpp_1\Downloads\2022-04-19_simu1_results.csv');  % ! heading problem
    % time = result(:,1)/1000;
    % v1 = result(:,2);
    % v1 = result(:,2);
    stimulus = true;
    try 
        pre_stim_ms =  750;
        post_stim_ms = 750; % ?
        timepoints_around_stimulus = (-pre_stim_ms:post_stim_ms);
        number_timepoints_around_stimulus = length(timepoints_around_stimulus);
        poststim_point = 250;  % ?
        index_poststim_point = pre_stim_ms+poststim_point+1;%+1 to skip the stimulation time

        %% Extract stimulus time series
        if stimulus
            tempo_data_file_path = fullfile(my_path,[simulation_name,stim_sizes_name{i},'_tempo.csv']);  
            stimulus_timeseries = readmatrix(tempo_data_file_path); %idx that stimuli are given
            stimulus_timeseries = floor(stimulus_timeseries);

            %turn from 0-based index in python to 1-based index in
            %matlab %! ????
            stimulus_timeseries = stimulus_timeseries + 1; % To align the index for Matlab style
        else
            %stimulus size is 0, then generate random series of
            %stimuli, and compute the plfs with respect to these
            %times
            no_stims = floor(length(raw_network_signal)/1000)-1;
            stimulusPeriodicity = 1000;
            jitter = 500;
            stimulus_timeseries = (1:no_stims)*stimulusPeriodicity - jitter/2 + unidrnd(jitter,1,no_stims);
        end

        %take out stimuli that are too close to the beginning / end
        %of the simulation, to allow for computation across the
        %entire prestimulus/poststimulus interval
        stimulus_timeseries = stimulus_timeseries(stimulus_timeseries>=pre_stim_ms & stimulus_timeseries<=length(raw_signal_to_process)-post_stim_ms);

        no_stims = length(stimulus_timeseries);  % number of stimulus

        %% Filter signal and extract phase and amplitude
        filtered = filter_fir(raw_signal_to_process,delta_hp,delta_lp,1000,2/delta_hp);

        filtered = filtered - mean(filtered); %!
        signal_phase = angle(hilbert(filtered)); % Get the phase for the whole time series
        signal_amplitude = abs(hilbert(filtered)); % filter the signal with a band

        %% Split phase and amplitudes timeseries by trial
        phase_per_trial = zeros(pre_stim_ms+post_stim_ms+1, no_stims);
        amplitude_per_trial  = zeros(pre_stim_ms+post_stim_ms+1, no_stims);

        for stim_idx = 1:no_stims
            stimulation_time = floor(stimulus_timeseries(stim_idx));
            phase_per_trial(:,stim_idx)  = signal_phase(stimulation_time-pre_stim_ms:(stimulation_time + post_stim_ms))';
            amplitude_per_trial(:,stim_idx)  = signal_amplitude(stimulation_time-pre_stim_ms:(stimulation_time + post_stim_ms))';
        end

        %% Compute pre-stimulus amplitude regulation for all timepoints around stimulus
        amplitude_regulation=zeros(size(timepoints_around_stimulus));

        prestimulus_amplitude_interval = [-150,-50];

        %calculate, for all trials, the amplitude in the prestimulus interval
        amplitude_prestim = zeros(size(stimulus_timeseries));
        for stim_idx=1:length(stimulus_timeseries)
            amplitude_prestim(stim_idx)=mean(signal_amplitude(stimulus_timeseries(stim_idx)+prestimulus_amplitude_interval(1):stimulus_timeseries(stim_idx)+prestimulus_amplitude_interval(2)));
        end

        %sort by prestimulus amplitude
        [~,stim_index] = sort(amplitude_prestim);

        %split into 10 percentiles, by prestimulus amplitude
        percentile_plf = zeros(10,1);                
        percentile_length = length(stim_index)/10;
        percentile_start_indexes = round(percentile_length * (0:9))+1;
        percentile_end_indexes = [percentile_start_indexes(2:10)-1 length(stim_index)];

        %compute prestimulus regulation for all timepoints around stimulus
        for t=1:length(timepoints_around_stimulus)
            %compute plf for trials in each bin
            for percentile_idx=1:10
                percentile_plf(percentile_idx)=circ_r(phase_per_trial(t,stim_index(percentile_start_indexes(percentile_idx):percentile_end_indexes(percentile_idx)))');
            end
            amplitude_regulation(t)= corr((1:10)',percentile_plf,'Type','Spearman');

            %also save the actual plfs (Phase Locking Factors) for the desired prestimulus
            %point
            if timepoints_around_stimulus(t)==poststim_point
                amplitude_percentile_plfs = percentile_plf;
            end
        end

        %% Compute prestimulus phase dependence for all timepoints around stimulus

        %make phase bins spaced at pi/16
        phase_regulation = zeros(1,number_timepoints_around_stimulus);
        phases_sort_time = -5;
        index_sort_time = pre_stim_ms+phases_sort_time;
        bins_phases = -pi:pi/16:pi;
        % each bin is assigned trials wih prestimulus phases within
        % a binwidth, centered on the bin
        bin_width = pi/4;
        for bin_idx = 1:length(bins_phases)
            orig = abs(circ_dist(phase_per_trial(index_sort_time,:),bins_phases(bin_idx)));
            pcks = orig < bin_width/2;
            no_pcks(bin_idx) = nnz(pcks);
        end
        min_pcks = min(no_pcks);

        phasebin_plfs_all = zeros(1,number_timepoints_around_stimulus);
        for bin_idx = 1:length(bins_phases)
            orig = abs(circ_dist(phase_per_trial(index_sort_time,:),bins_phases(bin_idx)));
            pcks = orig < bin_width/2;
            pcks = find(pcks);
            pcks = pcks(1:min_pcks);
            for stimulation_time = 1:number_timepoints_around_stimulus
                phasebin_plfs_all(bin_idx,stimulation_time) = circ_r(phase_per_trial(stimulation_time,pcks)');
            end
        end

        %% Compute phase regulation and plf for all timepoints around stimulus
        plfs = zeros(1,number_timepoints_around_stimulus);
        for stimulation_time = 1:number_timepoints_around_stimulus
            phase_regulation(stimulation_time) = circ_r(bins_phases',phasebin_plfs_all(:,stimulation_time));
            plfs(stimulation_time)=circ_r(phase_per_trial(stimulation_time,:)');
        end

        phasebin_plfs = phasebin_plfs_all(:,index_poststim_point);%contains the plfs for all phase bins at the poststimulus point
        ampdep_poststim = amplitude_regulation(index_poststim_point);
        phasedep_poststim = phase_regulation(index_poststim_point);
        plf_poststim = circ_r(phase_per_trial(index_poststim_point,:)');
        
        amplitude_regulations{i}=amplitude_regulation;
        phase_regulations{i}=phase_regulation;
        plfss{i}=plfs;
        amplitude_percentile_plfss{i}=amplitude_percentile_plfs;
        bins_phasess{i}=bins_phases;
        phasebin_plfss{i}=phasebin_plfs;
        PLF=[PLF;plf_poststim];
        %% save prestimulus regulation statistics
    %     save(fullfile(my_path,'prestim_dependence.mat'),'amplitude_regulation','phase_regulation','plfs','amplitude_percentile_plfs','bins_phases','phasebin_plfs');
    catch E
        ampdep_poststim=NaN;
        phasedep_poststim=NaN;
        plf_poststim=NaN;
        disp(getReport(E));
    end
end
T=table(Run_Number,amplitude_regulations,phase_regulations,plfss,amplitude_percentile_plfss,...
    bins_phasess,phasebin_plfss,Stimulus_size,PLF,...
    'RowNames', stim_sizes_name);
save(fullfile(my_path,[simulation_name,'_plf.mat']),'T');