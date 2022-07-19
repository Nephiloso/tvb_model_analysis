clear all;

% The script aims to compute the dynamic range to reproduce 1E,2B,3B using the data from TVB model.
% The script still need to be tested since the data is not ready yet.
% Most of the values of dynamic ranges are unrealistic => debug!

M= readtable('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\report\figures\critical_grid_analysis.csv');

addpath('support_functions_matlab');

data_path = 'D:\Downloads\data_intern_vu';
figures_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\report\figures\heatmaps';

model_analysis_code_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\Arthur 2020 code\model code\data_analysis';

addpath(genpath(model_analysis_code_path));

c_ee = 6:4:23;
c_ei = 6:4:23;
stimulus_size = [1; 10; 40; 100; 120; 150; 170];
stim_sizes_name = {'1'; '10'; '40'; '100'; '120'; '150'; '170'};
[X,Y] = meshgrid(c_ei,c_ee);

Z = zeros(length(c_ei),length(c_ee)); % dynamic range
V = zeros(length(c_ei),length(c_ee),length(stimulus_size));  % amp reg
U = zeros(length(c_ei),length(c_ee),length(stimulus_size));  % phs reg

file_date= '2022-07-13';
file_date1= '2022-07-14';
file_date2= '2022-07-15';
for i=c_ee
    for j=c_ei
        plf_trial=[];
        for k=1:length(stimulus_size)
            try
                stimulation_name = [file_date,'_c_ee',num2str(i),'_c_ei',num2str(j),'_stim_size',num2str(stimulus_size(k)),'_stim_results.csv'];
                tempo_name = [file_date,'_c_ee',num2str(i),'_c_ei',num2str(j),'_stim_size',num2str(stimulus_size(k)),'_tempo.csv'];
                simulation_data_file_path = fullfile(data_path,stimulation_name);
                tempo_data_file_path = fullfile(data_path,tempo_name);

                stimulus_timeseries = readmatrix(tempo_data_file_path); %idx that stimuli are given
                
            catch E
                try
                    stimulation_name = [file_date1,'_c_ee',num2str(i),'_c_ei',num2str(j),'_stim_size',num2str(stimulus_size(k)),'_stim_results.csv'];
                    tempo_name = [file_date1,'_c_ee',num2str(i),'_c_ei',num2str(j),'_stim_size',num2str(stimulus_size(k)),'_tempo.csv'];
                    simulation_data_file_path = fullfile(data_path,stimulation_name);
                    tempo_data_file_path = fullfile(data_path,tempo_name);

                    stimulus_timeseries = readmatrix(tempo_data_file_path); %idx that stimuli are given
                catch E
                    stimulation_name = [file_date2,'_c_ee',num2str(i),'_c_ei',num2str(j),'_stim_size',num2str(stimulus_size(k)),'_stim_results.csv'];
                    tempo_name = [file_date2,'_c_ee',num2str(i),'_c_ei',num2str(j),'_stim_size',num2str(stimulus_size(k)),'_tempo.csv'];
                    simulation_data_file_path = fullfile(data_path,stimulation_name);
                    tempo_data_file_path = fullfile(data_path,tempo_name);

                    stimulus_timeseries = readmatrix(tempo_data_file_path); %idx that stimuli are given
                end
            end
            
            stimulus_timeseries = floor(stimulus_timeseries);

            result = readmatrix(simulation_data_file_path);
            v1 = result(:,2);
            raw_signal_to_process = v1;
            plf_poststim=get_plf(raw_signal_to_process,stimulus_timeseries);
            plf_trial = [plf_trial,plf_poststim];
            
            phasedep_poststim=get_phasedep(raw_signal_to_process,stimulus_timeseries);
            
            idx = find(M.config_c_ee==i&M.config_c_ei==j);
            V(c_ee==i,c_ei==j,k) = plf_poststim;  % amp regulation
            U(c_ee==i,c_ei==j,k) = phasedep_poststim;  % phs regulation
            clear plf_poststim
        end
        
        crt_data=plf_trial(plf_trial>0);
        rn=plf_range(crt_data,stimulus_size);
        Z(c_ee==i,c_ei==j) = rn;
        M.dynamic_range(idx)=rn;
        
        clear plf_trial
        clear idx
    end
end
%% Plot the PLF dynamic range
fig = figure('color','w');
yvalues = cellstr (string(uint8(c_ee)));
xvalues = cellstr (string(uint8(c_ei)));
hm = heatmap(xvalues,yvalues,Z);
xlabel('W[I->E]')
ylabel('W[E->E]')
title('Dynamic range');
hm.YDisplayData = flipud(hm.YDisplayData);
colormap('jet');
caxis([0 10]);
colorbar

% saveas(fig,fullfile(figures_path,'fig_1E.fig'),'fig');

%% Plot the amplitude regulation
for k=length(stimulus_size)
    figure('color','w');
    yvalues = cellstr (string(uint8(c_ee)));
    xvalues = cellstr (string(uint8(c_ei)));
    hm = heatmap(xvalues,yvalues,V);
    xlabel('W[I->E]')
    ylabel('W[E->E]')
    title('PLF');
    hm.YDisplayData = flipud(hm.YDisplayData);
    colormap('jet');
    caxis([0 1]);
    title(['PLF - stim size: ',um2str(stimulus_size(k))]);
    colorbar
    saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_2B.fig']),'fig');
end
%% Plot the phase regulation
for u=length(stimulus_size)
    stimulus_size_to_analyze=stimulus_size(u);
    figure('color','w');
    yvalues = cellstr (string(uint8(c_ee)));
    xvalues = cellstr (string(uint8(c_ei)));
    hm = heatmap(xvalues,yvalues,U(:,:,stimulus_size_to_analyze));
    xlabel('W[I->E]')
    ylabel('W[E->E]')
    title('PLF');
    hm.YDisplayData = flipud(hm.YDisplayData);
    colormap('jet');
    caxis([0 1]);
    colorbar
    saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_3B.fig']),'fig');
end

%% save the data
M1.config_c_ee = M.config_c_ee;
M1.config_c_ei = M.config_c_ei;
M1.dfa_theta = M.dfa_all_theta;
M1.dynamic_range = M.dynamic_range;
% save(fullfile(figures_path,'critical_grid_plf_analysis.m'),'M1');
% save(fullfile(figures_path,'1E_data.m'),'Z');
% save(fullfile(figures_path,'2B_data.m'),'V');
% save(fullfile(figures_path,'3B_data.m'),'U');
function plf_poststim=get_plf(raw_signal_to_process,stimulus_timeseries)
    delta_lp = 4;
    delta_hp = 0.5;
    try 
        pre_stim_ms =  750;
        post_stim_ms = 750;
        poststim_point = 250;
        index_poststim_point = pre_stim_ms+poststim_point+1;%+1 to skip the stimulation time
        Fs=1000;
        stimulus_timeseries = stimulus_timeseries(stimulus_timeseries>=pre_stim_ms & stimulus_timeseries<=length(raw_signal_to_process)-post_stim_ms);

        no_stims = length(stimulus_timeseries);  % number of stimulus

        %% Filter signal and extract phase and amplitude
        filtered = filter_fir(raw_signal_to_process,delta_hp,delta_lp,Fs,2/delta_hp);

        filtered = filtered - mean(filtered);
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
                
        plf_poststim = circ_r(phase_per_trial(index_poststim_point,:)');
    catch E
        plf_poststim=NaN;
        disp(getReport(E));
    end
end

function rn=plf_range(crt_data,crt_stims)
    fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
    fsigm_inv = @(param,xval) param(3)-reallog((param(2)-param(1))./(xval-param(1)) - 1) * 1/param(4);

    try
        params = sigm_fit(reallog(crt_stims),crt_data);

        x = reallog(crt_stims);
        fitDats = fsigm(params,x);

        mnFit = min(fitDats);
        mxFit = max(fitDats);
        tenPerc = (mxFit - mnFit)*0.1;
        ninetyPerc = mxFit - tenPerc;
        tenPerc = tenPerc + mnFit;

        xinv_10perc = fsigm_inv(params,tenPerc);
        xinv_90perc = fsigm_inv(params,ninetyPerc);

        %this ends up empty if the 
        if isempty(xinv_10perc) | isempty(xinv_90perc)
            rn = NaN;
        else
            if xinv_10perc > xinv_90perc
                rn = NaN;
            else
                rn = xinv_90perc-xinv_10perc;
            end
        end

    catch E
        disp(getReport(E));
        rn = NaN;
    end
end

function phasedep_poststim=get_phasedep(raw_signal_to_process,stimulus_timeseries)
    delta_lp = 4;
    delta_hp = 0.5;
    pre_stim_ms =  750;
    post_stim_ms = 750; % ?
    timepoints_around_stimulus = (-pre_stim_ms:post_stim_ms);
    number_timepoints_around_stimulus = length(timepoints_around_stimulus);
    poststim_point = 250;  % ?
    index_poststim_point = pre_stim_ms+poststim_point+1;
    
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
    for stimulation_time = 1:number_timepoints_around_stimulus
        phase_regulation(stimulation_time) = circ_r(bins_phases',phasebin_plfs_all(:,stimulation_time));
    end
    phasedep_poststim = phase_regulation(index_poststim_point);
end