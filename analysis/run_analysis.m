%     Originally created by Arthur-Ervin Avramiea (2020), arthur.avramiea@gmail.com

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function run_analysis( my_path )
    try
        if my_path(length(my_path))~='/' && my_path(length(my_path))~='\'
            my_path = [my_path filesep];
        end
        
        disp(['analysing signal at path' my_path]);
        
 		parameters_dict=read_dictionary_from_file(strcat(my_path,'Simulation.Parameters'));
       
        %% load analysis parameters
 		
        dfa_calc_smallest_window = parameters_dict('DFA calc smallest window');
 		dfa_calc_largest_window = parameters_dict('DFA calc largest window');
 		dfa_calc_window_overlap = parameters_dict('DFA calc window overlap');
 		dfa_fit_smallest_window = parameters_dict('DFA fit smallest window');
 		dfa_fit_largest_window = parameters_dict('DFA fit largest window');
 
        delta_hp = parameters_dict('delta_hp');
        delta_lp = parameters_dict('delta_lp');
        theta_hp = parameters_dict('theta_hp');
        theta_lp = parameters_dict('theta_lp');
        alpha_hp = parameters_dict('alpha_hp');
        alpha_lp = parameters_dict('alpha_lp');
        beta_hp = parameters_dict('beta_hp');
        beta_lp = parameters_dict('beta_lp');
        gamma_hp = parameters_dict('gamma_hp');
        gamma_lp = parameters_dict('gamma_lp');
        
 		kappa_avalanche_size_min = parameters_dict('Kappa avalanche size min');
 		kappa_avalanche_size_max = parameters_dict('Kappa avalanche size max');
 		kappa_bins = parameters_dict('Kappa bins');
        kappa_threshold_type = parameters_dict('Kappa threshold type');
        
        external_input=parameters_dict('External input');
        stimulus_size=parameters_dict('Stimulus size');
        
        add_whitenoise_first=parameters_dict('Add white noise');
        white_noise_sd=parameters_dict('White noise sd');
        
        Fs = 1000;

        simulation_data_file_path = fullfile(my_path,'SimulationData.h5');  
        excitatory_spikes_time_series = h5read(simulation_data_file_path,'/data/main_network_excitatory_spikes_timeseries');
        inhibitory_spikes_time_series = h5read(simulation_data_file_path,'/data/main_network_inhibitory_spikes_timeseries');
        raw_network_signal = excitatory_spikes_time_series + inhibitory_spikes_time_series;
		
        if isempty(raw_network_signal)
            return;
        end
        
        %shuffle random number generator
		rng('shuffle');
        
        if add_whitenoise_first
            noise = randn(size(raw_network_signal,1),1)*white_noise_sd;
            raw_signal_to_process = raw_network_signal + noise;
        else
            raw_signal_to_process = raw_network_signal;
        end
        
        %% get amplitude envelope in different bands
        AmplitudeEnvelope_alpha = get_amplitude_envelope(raw_signal_to_process, Fs, alpha_hp, alpha_lp);
        AmplitudeEnvelope_beta  = get_amplitude_envelope(raw_signal_to_process, Fs, beta_hp, beta_lp);
        AmplitudeEnvelope_gamma = get_amplitude_envelope(raw_signal_to_process, Fs, gamma_hp, gamma_lp);
        AmplitudeEnvelope_delta = get_amplitude_envelope(raw_signal_to_process, Fs, delta_hp, delta_lp);
        AmplitudeEnvelope_theta = get_amplitude_envelope(raw_signal_to_process, Fs, theta_hp, theta_lp); 
        
        %% compute activity (average number of spikes per timestep)
        mean_act = mean(raw_network_signal);
        
        %% DFA calculation
        try
            [x_alpha,y_alpha]=calculate_DFA_points(AmplitudeEnvelope_alpha,Fs,[dfa_calc_smallest_window,dfa_calc_largest_window],dfa_calc_window_overlap);
            DFA_exp_alpha=fit_DFA(x_alpha',y_alpha,[dfa_fit_smallest_window,dfa_fit_largest_window],Fs);

            [x_beta,y_beta]=calculate_DFA_points(AmplitudeEnvelope_beta,Fs,[dfa_calc_smallest_window,dfa_calc_largest_window],dfa_calc_window_overlap);
            DFA_exp_beta=fit_DFA(x_beta',y_beta,[dfa_fit_smallest_window,dfa_fit_largest_window],Fs);

            [x_delta,y_delta]=calculate_DFA_points(AmplitudeEnvelope_delta,Fs,[dfa_calc_smallest_window,dfa_calc_largest_window],dfa_calc_window_overlap);
            DFA_exp_delta=fit_DFA(x_delta',y_delta,[dfa_fit_smallest_window,dfa_fit_largest_window],Fs);

            [x_theta,y_theta]=calculate_DFA_points(AmplitudeEnvelope_theta,Fs,[dfa_calc_smallest_window,dfa_calc_largest_window],dfa_calc_window_overlap);
            DFA_exp_theta=fit_DFA(x_theta',y_theta,[dfa_fit_smallest_window,dfa_fit_largest_window],Fs);

            [x_gamma,y_gamma]=calculate_DFA_points(AmplitudeEnvelope_gamma,Fs,[dfa_calc_smallest_window,dfa_calc_largest_window],dfa_calc_window_overlap);
            DFA_exp_gamma=fit_DFA(x_gamma',y_gamma,[dfa_fit_smallest_window,dfa_fit_largest_window],Fs);
        catch E
            DFA_exp_alpha=NaN;
            disp(getReport(E));
        end
        
        %% Kappa calculation
        try
            % avalanche analysis is performed on the raw spike time
            % series, without any noise added
            if kappa_threshold_type==1
                threshold=median(raw_network_signal)/2;
            else
                threshold=0;
            end
            [Kappa_size,Kappa_dur,~] = calculate_kappa(raw_network_signal',threshold,kappa_bins,kappa_avalanche_size_min,kappa_avalanche_size_max);
        catch E
            Kappa_size=NaN;
            Kappa_dur=NaN;
            disp(getReport(E));
        end
        
        %% Compute power and amplitude
        nfft=2^11; %number of fast fourier transforms, the higher this number, the more the frequency resolution of the spectrum goes up

        % define intervals at which power is calculated
        interval_Hz(1,:)=[delta_hp delta_lp];
        interval_Hz(2,:)=[theta_hp theta_lp];
        interval_Hz(3,:)=[alpha_hp alpha_lp];
        interval_Hz(4,:)=[beta_hp beta_lp];
        interval_Hz(5,:)=[gamma_hp gamma_lp];

        min_freq = min(interval_Hz(:));
        max_freq = max(interval_Hz(:));

        % calculated integrated and normalized power, and amplitude for all channels
        [power,f]=pwelch(raw_signal_to_process,hamming(nfft),[],[],Fs);
        amplitude=sqrt(power); % transform power to amplitude

        for interval_idx=1:size(interval_Hz,1)
            interval{interval_idx}=find(f>interval_Hz(interval_idx,1)&f<interval_Hz(interval_idx,2));
        end

        nr_interv=length(interval);

        absolute_amplitude=zeros(nr_interv,1);
        relative_amplitude=zeros(nr_interv,1);
        for j=1:nr_interv
            absolute_amplitude(j,1) = mean(amplitude(interval{j}));
            relative_amplitude(j,1) = sum(amplitude(interval{j}))/sum ( amplitude ( f>min_freq & f<max_freq )) ;
        end

        absolute_power=zeros(nr_interv,1);
        relative_power=zeros(nr_interv,1);
        for j=1:nr_interv
            absolute_power(j,1) = mean(power(interval{j}));
            relative_power(j,1) = sum(power(interval{j}))/sum( power( f>min_freq & f<max_freq ) );
        end
        
        
        %% Run stimulus response analysis
        if external_input
            try 
                pre_stim_ms =  750;
                post_stim_ms = 750;
                timepoints_around_stimulus = (-pre_stim_ms:post_stim_ms);
                number_timepoints_around_stimulus = length(timepoints_around_stimulus);
                poststim_point = 150;
                index_poststim_point = pre_stim_ms+poststim_point+1;%+1 to skip the stimulation time
                
                %% Extract stimulus time series
                if stimulus_size > 0
                    stimulus_timeseries = h5read(simulation_data_file_path,'/data/external_input_spike_times');
                    stimulus_timeseries = floor(stimulus_timeseries);
                    
                    %turn from 0-based index in python to 1-based index in
                    %matlab
                    stimulus_timeseries = stimulus_timeseries + 1;
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

                no_stims = length(stimulus_timeseries);

                %% Filter signal and extract phase and amplitude
                filtered = filter_fir(raw_signal_to_process,alpha_hp,alpha_lp,1000,2/alpha_hp);

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
                
                %% Compute pre-stimulus amplitude regulation for all timepoints around stimulus
                amplitude_regulation=zeros(size(timepoints_around_stimulus));
                
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
                percentile_plf = zeros(10,1);                
                percentile_length = length(stim_index)/10;
                percentile_start_indexes = round(percentile_length * (0:9))+1;
                percentile_end_indexes = [percentile_start_indexes(2:10)-1 length(stim_index)];
                
                %compute prestimulus regulation for all timepoints around
                %stimulus
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
                
                %% save prestimulus regulation statistics
                save(strcat(my_path,'prestim_dependence.mat'),'amplitude_regulation','phase_regulation','plfs','amplitude_percentile_plfs','bins_phases','phasebin_plfs');
            catch E
                ampdep_poststim=NaN;
                phasedep_poststim=NaN;
                plf_poststim=NaN;
                disp(getReport(E));
            end
        end
        
        %% Write out calculated statistics
        activity_statistics_file_path = strcat(my_path,'ActivityStatistics.txt') ;
        
        if isfile( activity_statistics_file_path )
            already_existing_parameters = read_dictionary_from_file(activity_statistics_file_path);
        else
            already_existing_parameters = containers.Map();
        end
        
        new_params_to_write = containers.Map();
		
        new_params_to_write('DFA_alpha')=DFA_exp_alpha;
        new_params_to_write('DFA_beta')=DFA_exp_beta;
        new_params_to_write('DFA_delta')=DFA_exp_delta;
        new_params_to_write('DFA_gamma')=DFA_exp_gamma;
        new_params_to_write('DFA_theta')=DFA_exp_theta;

        new_params_to_write('delta_absolute_power')=absolute_power(1,1);
        new_params_to_write('theta_absolute_power')=absolute_power(2,1);
        new_params_to_write('alpha_absolute_power')=absolute_power(3,1);
        new_params_to_write('beta_absolute_power')=absolute_power(4,1);
        new_params_to_write('gamma_absolute_power')=absolute_power(5,1);

        new_params_to_write('delta_relative_power')=relative_power(1,1);
        new_params_to_write('theta_relative_power')=relative_power(2,1);
        new_params_to_write('alpha_relative_power')=relative_power(3,1);
        new_params_to_write('beta_relative_power')=relative_power(4,1);
        new_params_to_write('gamma_relative_power')=relative_power(5,1);

        new_params_to_write('delta_absolute_amplitude')=absolute_amplitude(1,1);
        new_params_to_write('theta_absolute_amplitude')=absolute_amplitude(2,1);
        new_params_to_write('alpha_absolute_amplitude')=absolute_amplitude(3,1);
        new_params_to_write('beta_absolute_amplitude')=absolute_amplitude(4,1);
        new_params_to_write('gamma_absolute_amplitude')=absolute_amplitude(5,1);

        new_params_to_write('delta_relative_amplitude')=relative_amplitude(1,1);
        new_params_to_write('theta_relative_amplitude')=relative_amplitude(2,1);
        new_params_to_write('alpha_relative_amplitude')=relative_amplitude(3,1);
        new_params_to_write('beta_relative_amplitude')=relative_amplitude(4,1);
        new_params_to_write('gamma_relative_amplitude')=relative_amplitude(5,1);

        new_params_to_write('Kappa_size')=Kappa_size;
        new_params_to_write('Kappa_dur')=Kappa_dur;
        
        if external_input
            new_params_to_write('PLF')=plf_poststim;
            new_params_to_write('Prestim_ampdep')=ampdep_poststim;
            new_params_to_write('Prestim_phasedep')=phasedep_poststim;
        end
        
        new_params_to_write('mean_act')=mean_act;
        
        %loaded parameters from existing file
        params_to_write = already_existing_parameters;
        %overwrite with newly calculated
        keys_new_params_to_write = keys(new_params_to_write);
        for i=1:length(keys_new_params_to_write)
            crt_key_to_write = keys_new_params_to_write{i};
            params_to_write(crt_key_to_write)=new_params_to_write(crt_key_to_write);
        end
        
        write_params(my_path,'ActivityStatistics.txt',params_to_write);
    catch E
        disp(getReport(E));
    end
end
