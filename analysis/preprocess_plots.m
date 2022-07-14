% The script aims to calculate the dfa and kappa for the signal created by the virtual brain.
% The DFA analysis has been moved to Python. This script is abandoned and saved as a backup.

parameters_dict=read_dictionary_from_file('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\V1criticality\analysis\Simulation.Parameters');

% add the analysis folder to path
addpath('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\Arthur 2020 code\model code\data_analysis')

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

calc_kappa=parameters_dict('Calculate Kappa');
kappa_avalanche_size_min = parameters_dict('Kappa avalanche size min');
kappa_avalanche_size_max = parameters_dict('Kappa avalanche size max');
kappa_bins = parameters_dict('Kappa bins');
kappa_threshold_type = parameters_dict('Kappa threshold type');

calc_DFA=parameters_dict('Calculate DFA');
calc_amplitude = parameters_dict('Calculate amplitude');
add_whitenoise_first=parameters_dict('Add white noise');
white_noise_sd=parameters_dict('White noise sd');

random_seed = 6;

Fs = 1000;  % monitors-period: ms => Fs: Hz

result = readmatrix('C:\Users\wpp_1\Downloads\2022-06-27_c_ee14_c_ei14_results.csv');  % ! heading problem

time = result(:,1)/1000;
v1 = result(:,2);
% v2 = result(:,3);

% raw_signal_to_process = table2array(raw_signal_to_process);
white_noise_sd = 0.1;
if add_whitenoise_first
    noise = randn(size(v1,1),1)*white_noise_sd;
    raw_signal_to_process = v1 + noise;
else
    raw_signal_to_process = v1;
end

%% Visualization of the raw data
figure('Name','signal profile - V1','Position',[50,50,1020,500])
plot(time(15000:25000),v1(15000:25000),'color',[0.18 0.5 0.1]);
title('Raw signal 5s')
xlabel('time/s')
ylabel('Average firing rate')

figure('Name','signal profile - V1','Position',[50,50,1120,640])
plot(time(15000:55000),v1(15000:55000));
title('Raw signal 40s - spatial average')
xlabel('time/s')
ylabel('populational firing rate')

figure('Name','signal +noise - V1','Position',[50,50,1120,640])
plot(time(15000:55000),raw_signal_to_process(15000:55000));
title('Raw signal 40s - spatial average')
xlabel('time/s')
ylabel('populational firing rate')
%% Power spectrum
figure()
nfft            = 8000 ; % for freq resolution = 0.125Hz (fs(1000)/0.125 = 8000) \ Controls freq resolution of psd (most commonly preferred resolution is 0.1Hz)
% raw data to be in Nx1 format; hamming windowed PSD (other methods also available)
[pxx,f]         = pwelch(v1,hamming(nfft),0,nfft,Fs,'psd');

% plot
plot(f,10*log10(pxx),'g','linewidth',1.8,'color',[0.18,0.5,0.1])
xlim([1 40]) % adjust for visualization
legend off
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
title('Power spectral density')
%% get amplitude envelope in different bands  input?spike series

AmplitudeEnvelope_delta = get_amplitude_envelope(raw_signal_to_process, Fs, delta_hp, delta_lp);

delta_filtered = filter_fir(raw_signal_to_process(:,:),delta_hp,delta_lp,Fs,2/delta_hp);
figure('Name','signal profile - V1','Position',[50,50,1020,250])
plot(time(15000:28000),delta_filtered(15000:28000),'linewidth',0.7,'color',[0 0.78 0.55]);
hold on
plot(time(15000:28000),AmplitudeEnvelope_delta(15000:28000),'linewidth',1.2,'color',[0.18 0.5 0.1]);
xlabel('time/s')
ylabel('Average firing rate')