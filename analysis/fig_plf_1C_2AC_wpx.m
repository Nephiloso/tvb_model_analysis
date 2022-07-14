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


%%% !!! Test: add white noise or not!
addpath('support_functions_matlab');

data_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\optimiz_stim';
figures_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\report\figures';
model_analysis_code_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\Arthur 2020 code\model code\data_analysis';

addpath(genpath(model_analysis_code_path));

noncrit_run_number	= 0;

crit_run_number = 1;

stimulus_size_to_analyze = 40;

non_cri_simulation_name = '2022-07-13_c_ee6_c_ei22_stim_size';
cri_simulation_name = '2022-07-05_c_ee11.57199543729056_stim_size';

%plot plf for critical networks at different stimulus sizes
plf_crit_0_stim = get_plf_all_trials(data_path,...
    [cri_simulation_name,'0','_stim_results.csv'],[cri_simulation_name,'0','_tempo.csv']);
plf_crit_10_stim = get_plf_all_trials(data_path,...
    [cri_simulation_name,'10','_stim_results.csv'],[cri_simulation_name,'10','_tempo.csv']);
plf_crit_100_stim = get_plf_all_trials(data_path,...
    [cri_simulation_name,'100','_stim_results.csv'],[cri_simulation_name,'100','_tempo.csv']);
plf_crit_170_stim = get_plf_all_trials(data_path,...
    [cri_simulation_name,'170','_stim_results.csv'],[cri_simulation_name,'170','_tempo.csv']);

fig = figure('color','w');
set(fig,'Position',[0 0 600 500]);
plot(-0.75:0.001:0.75,plf_crit_0_stim,'LineWidth',2,'Color',[0 0 0]);
hold on;
plot(-0.75:0.001:0.75,plf_crit_10_stim,'LineWidth',2,'Color',[0.0351 0.664 0.203]);
plot(-0.75:0.001:0.75,plf_crit_100_stim,'LineWidth',2,'Color',[0.05 0.945 0.289]);
plot(-0.75:0.001:0.75,plf_crit_170_stim,'LineWidth',2,'Color',[0.6 0.98 0.7]);
xlim([-0.75 0.75]);
ylim([0 1]);
legend({'0','10','100','170'},'Location','northwest');
xticks([-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
xticklabels({'-0.6','','','0','','','0.6'});
set(gca,'fontsize', 16);
xlabel({'Time since stimulation';'(seconds)'});
ylabel('PLF');
saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_1C.fig']),'fig');

%plot plf for subcrit/crit/supercrit at same stimulus size(5)
plf_noncrit_150_stim = get_plf_all_trials(data_path,...
    [non_cri_simulation_name,num2str(stimulus_size_to_analyze),'_stim_results.csv'],[non_cri_simulation_name,num2str(stimulus_size_to_analyze),'_tempo.csv']);
plf_crit_150_stim = get_plf_all_trials(data_path,...
    [cri_simulation_name,num2str(stimulus_size_to_analyze),'_stim_results.csv'],[cri_simulation_name,num2str(stimulus_size_to_analyze),'_tempo.csv']);
% plf_supercrit_5_stim = get_plf_all_trials(fullfile(data_path,'runs','without individual spikes',sprintf('EC%.2f.IC%.2f.Stim%i.Run%i',supercrit_exc_conn,supercrit_inh_conn,stimulus_size_to_analyze,supercrit_run_number)));

fig = figure('color','w');
set(fig,'Position',[0 0 600 500]);
plot(-0.75:0.001:0.75,plf_noncrit_150_stim,'LineWidth',2,'Color',[0 0.5 1]);
hold on;
plot(-0.75:0.001:0.75,plf_crit_150_stim,'LineWidth',2,'Color',[0.22 1 0.22]);
% plot(-0.75:0.001:0.75,plf_crit_150_stim,'LineWidth',1.5,'Color',[1 0.22 0.22]);
xlabel({'Time since stimulation';'(seconds)'});
ylabel('PLF');
ylim([0 0.4]);
yticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4]);
yticklabels({'0','','0.1','','0.2','','0.3','','0.4'});
xlim([-0.75 0.75]);
xticks([-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
xticklabels({'-0.6','','','0','','','0.6'});
set(gca,'fontsize', 16);
legend({'Noncritical','Critical'},'Location','northwest');
saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_2A.fig']),'fig');

%plot plf at different pre-stimulus amplitude percentiles for critical network at same stimulus
%size(150)
fig = plot_plf_percentiles(data_path,...
    [cri_simulation_name,num2str(stimulus_size_to_analyze),'_stim_results.csv'],[cri_simulation_name,num2str(stimulus_size_to_analyze),'_tempo.csv']);
saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_2C.fig']),'fig');

function fig = plot_plf_percentiles(folder_path, simu_file_name, stimu_file_name)

    [time_series,stimulus_timeseries] = get_timeseries(folder_path, simu_file_name, stimu_file_name);
%     time_series_noise = add_white_noise(time_series,3);
    time_series_noise = time_series;

    delta_lp = 4;
    delta_hp = 0.5;

    pre_stim_ms =  750;
    post_stim_ms = 750;
    timepoints_around_stimulus = (-pre_stim_ms:post_stim_ms);

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

    for interval_idx = 1:no_stims
        stimulation_time = floor(stimulus_timeseries(interval_idx));
        phase_per_trial(:,interval_idx)  = signal_phase(stimulation_time-pre_stim_ms:(stimulation_time + post_stim_ms))';
        amplitude_per_trial(:,interval_idx)  = signal_amplitude(stimulation_time-pre_stim_ms:(stimulation_time + post_stim_ms))';
    end

    %% Compute plf for different pre-stimulus amplitude percentiles

    prestimulus_amplitude_interval = [-150,-50];

    %calculate, for all trials, the amplitude in the
    %prestimulus interval
    amplitude_prestim = zeros(size(stimulus_timeseries));
    for interval_idx=1:length(stimulus_timeseries)
        amplitude_prestim(interval_idx)=mean(signal_amplitude(stimulus_timeseries(interval_idx)+prestimulus_amplitude_interval(1):stimulus_timeseries(interval_idx)+prestimulus_amplitude_interval(2)));
    end

    %sort by prestimulus amplitude
    [~,stim_index] = sort(amplitude_prestim);

    %split into 10 percentiles, by prestimulus amplitude
    percentile_plf = zeros(10,length(timepoints_around_stimulus));                
    percentile_length = length(stim_index)/10;
    percentile_start_indexes = round(percentile_length * (0:9))+1;
    percentile_end_indexes = [percentile_start_indexes(2:10)-1 length(stim_index)];

    %compute PLF across trials in all percentile bins, at all time points
    for t=1:length(timepoints_around_stimulus)
        %compute plf for trials in each bin
        for interval_idx=1:10
            percentile_plf(interval_idx,t)=circ_r(phase_per_trial(t,stim_index(percentile_start_indexes(interval_idx):percentile_end_indexes(interval_idx)))');
        end
    end

    fig = figure('color','w');
    set(fig,'Position',[0 0 600 500]);
    plot(-0.75:0.001:0.75,percentile_plf(1,:),'LineWidth',2,'Color',[0 0 0]);
    hold on;
    plot(-0.75:0.001:0.75,percentile_plf(5,:),'LineWidth',2,'Color',[0.0351 0.664 0.203]);
    plot(-0.75:0.001:0.75,percentile_plf(10,:),'LineWidth',2,'Color',[0.6 0.98 0.7]);

    xticks([-0.6 -0.4 -0.2 0 0.2 0.4 0.6]);
    xticklabels({'-0.6','','','0','','','0.6'});
    yticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1]);
    yticklabels({'0','','0.1','','0.2','','0.3','','0.4','','0.5','','0.6','','0.7','','0.8','','0.9','','1.0'});
    xlabel({'Time since stimulation','(seconds)'});
    ylabel('PLF');
    legend('0-10th percentile','40-50th percentile','90-100th percentile','Location','northwest');

    ylim([0 1]);
    xlim([-0.75 0.75]);

    set(gca,'fontsize', 16);
end


function plfs = get_plf_all_trials(folder_path, simu_file_name, stimu_file_name)

    [time_series,stimulus_timeseries] = get_timeseries(folder_path, simu_file_name, stimu_file_name);
%     time_series_noise = add_white_noise(time_series,3);
    time_series_noise = time_series;
    delta_lp = 4;
    delta_hp = 0.5;

    pre_stim_ms =  750;
    post_stim_ms = 750;
    timepoints_around_stimulus = (-pre_stim_ms:post_stim_ms);
    number_timepoints_around_stimulus = length(timepoints_around_stimulus);
    
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

    plfs = zeros(1,number_timepoints_around_stimulus);
    for stimulation_time = 1:number_timepoints_around_stimulus
        plfs(stimulation_time)=circ_r(phase_per_trial(stimulation_time,:)');
    end
end