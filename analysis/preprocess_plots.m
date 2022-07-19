% The script aims to calculate the dfa and kappa for the signal created by the virtual brain.
% The DFA analysis has been moved to Python. This script is abandoned and saved as a backup.
clear all
parameters_dict=read_dictionary_from_file('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\V1criticality\analysis\Simulation.Parameters');

% add the analysis folder to path
addpath('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\Arthur 2020 code\model code\data_analysis')

figures_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\report\figures\preprocess';

data_folder = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\preprocess';
data_prefix = '2022-06-27_c_ee6_c_ei22';
% 2022-06-27_c_ee6_c_ei22_results.csv
% 2022-07_15_randm2_results.csv
% 2022-06-26_c_ee11.57199543729056_c_ei10.307919639315172_dt1_results.csv
data_name = [data_prefix,'_results.csv'];
filtered_name = [data_prefix,'_filtered0.5_4.csv'];
hilbert_name = [data_prefix,'_hilbert0.5_4.csv'];
psd_name = [data_prefix,'_psd.csv'];

delta_lp = 4;
delta_hp = 0.5;

Fs = 1000;  % monitors-period: ms => Fs: Hz

deep = [[0.18,0.5,0.1];[0.1,0.1,0.44];[0 0 0]];
shallow = [[0 0.78 0.55];[0 0.5 1];[0.5 0.5 0.5]];
mode = 2;
%% Visualization of the raw data
result = readmatrix(fullfile(data_folder,data_name));  % ! heading problem
time = result(:,1)/1000; % time in sec
v1 = result(:,2);

% raw_signal_to_process = table2array(raw_signal_to_process);
white_noise_sd = 0.1;
add_whitenoise_first = false;
if add_whitenoise_first
    noise = randn(size(v1,1),1)*white_noise_sd;
    raw_signal_to_process = v1 + noise;
else
    raw_signal_to_process = v1;
end

starting_pt = 15000; % data point
simu_length = 3; % second
end_pt = starting_pt+simu_length*Fs;
fig = figure('Name','signal profile - V1','Position',[50,50,1020,250],'color','w');
plot(time(starting_pt:end_pt),v1(starting_pt:end_pt),'color',shallow(mode,:));%[0.18 0.5 0.1]);
%title('Raw signal 5s')
% ylim([-1 1])
xlabel('time/s')
ylabel('Average firing rate')
% saveas(fig,fullfile(figures_path,[data_prefix,'raw.fig']),'fig');
% saveas(fig,fullfile(figures_path,[data_prefix,'raw.svg']),'svg');
% figure('Name','signal profile - V1','Position',[50,50,1120,640])
% plot(time(15000:55000),v1(15000:55000));
% %title('Raw signal 40s - spatial average')
% xlabel('time/s')
% ylabel('populational firing rate')
%% get amplitude envelope in different bands
starting_pt=15000;
simu_length=3;
end_pt = starting_pt+simu_length*Fs;
filtered_result = readmatrix(fullfile(data_folder,filtered_name));
delta_filtered = filtered_result(:,2);
Envelope_delta_data = readmatrix(fullfile(data_folder,hilbert_name));
AmplitudeEnvelope_delta = Envelope_delta_data(:,2);

figure('Name','filtered signal in delta band','Position',[50,50,1020,250],'color','w')
plot(time(starting_pt:end_pt),delta_filtered(starting_pt:end_pt),'linewidth',0.7,'color',shallow(mode,:));%[0 0.78 0.55]);
hold on
plot(time(starting_pt:end_pt),AmplitudeEnvelope_delta(starting_pt:end_pt),'linewidth',1.2,'color',deep(mode,:));
xlabel('time/s')
ylabel('Average firing rate')
%%%%%%%%%%%%%%%%%%% Filtered in matlab %%%%%%%%%%%%%%%%%%%%
% cent_signal = raw_signal_to_process-mean(raw_signal_to_process);
% AmplitudeEnvelope_delta = get_amplitude_envelope(cent_signal, Fs, delta_hp, delta_lp);
% delta_filtered = filter_fir(cent_signal(:,:),delta_hp,delta_lp,Fs,2/delta_hp);
% fig=figure('Name','filtered signal in delta band - 3s','Position',[50,50,1020,250],'color','w');
% plot(time(starting_pt:end_pt),delta_filtered(starting_pt:end_pt),'linewidth',1,'color',shallow(mode,:));
% hold on
% plot(time(starting_pt:end_pt),AmplitudeEnvelope_delta(starting_pt:end_pt),'linewidth',1.2,'color',deep(mode,:));
% xlabel('time/s')
% ylabel('Average firing rate')
% saveas(fig,fullfile(figures_path,[data_prefix,'filteres_short.fig']),'fig');
% saveas(fig,fullfile(figures_path,[data_prefix,'filteres_short.svg']),'svg');
% starting_pt=15000;
% simu_length=50;
% end_pt = starting_pt+simu_length*Fs;
% delta_filtered = filter_fir(cent_signal(:,:),delta_hp,delta_lp,Fs,2/delta_hp);
% fig=figure('Name','filtered signal in delta band - 50s','Position',[50,50,1020,250],'color','w');
% plot(time(starting_pt:end_pt),delta_filtered(starting_pt:end_pt),'linewidth',1,'color',shallow(mode,:));
% hold on
% plot(time(starting_pt:end_pt),AmplitudeEnvelope_delta(starting_pt:end_pt),'linewidth',1.2,'color',deep(mode,:));
% xlabel('time/s')
% ylabel('Average firing rate')
% saveas(fig,fullfile(figures_path,[data_prefix,'filteres_long.fig']),'fig');
% saveas(fig,fullfile(figures_path,[data_prefix,'filteres_long.svg']),'svg');
%% Visualization of the Power spectrum
psd_result = readmatrix(fullfile(data_folder,psd_name));
f = psd_result(:,1);
pxx = psd_result(:,2);

% plot
fig =figure('Name','Power Spectrum Density','Position',[50,50,470,250],'color','w');
plot(f,10*log10(pxx),'g','linewidth',1.8,'color',deep(mode,:))
xlim([1 40]) % adjust for visualization
legend off
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
%title('Power spectral density')
% saveas(fig,fullfile(figures_path,[data_prefix,'psd.fig']),'fig');
% saveas(fig,fullfile(figures_path,[data_prefix,'psd.svg']),'svg');