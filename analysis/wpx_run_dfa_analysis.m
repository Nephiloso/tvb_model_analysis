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
data_prefix = '2022-06-26_c_ee11.57199543729056_c_ei10.307919639315172_dt1';
data_name = [data_prefix,'_results.csv'];
% result = readmatrix('C:\Users\wpp_1\Downloads\2022-06-27_c_ee14_c_ei14_results.csv');  % ! heading problem
data_folder = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\preprocess';
result = readmatrix(fullfile(data_folder,data_name));
time = result(:,1)/1000;
v1 = result(:,2);
% v2 = result(:,3);

% raw_signal_to_process = table2array(raw_signal_to_process);
% add_whitenoise_first = false;
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
% title('Raw signal 5s')
xlabel('time/s')
ylabel('Average firing rate')
% h3 = subplot(3,1,3);

figure('Name','signal profile - V1','Position',[50,50,1120,640])
% subplot(3,1,1)
plot(time(15000:55000),v1(15000:55000));
title('Raw signal 40s - spatial average')
xlabel('time/s')
ylabel('populational firing rate')

figure('Name','signal +noise - V1','Position',[50,50,1120,640])
% subplot(3,1,1)
plot(time(15000:55000),raw_signal_to_process(15000:55000));
title('Raw signal 40s - spatial average')
xlabel('time/s')
ylabel('populational firing rate')

% alpha_filtered = bandpass(raw_signal_to_process,[alpha_hp alpha_lp],Fs);
% % subplot(3,1,2)
% plot(h3, time(1:50000),alpha_filtered(1:50000));
% title(h3, 'filtered signal(8-16Hz) 25s - spatial average')
% xlabel(h3, 'time/s')

% 
% subplot(3,1,3)
% plot(time(15000:105000),v1(15000:105000));
% title('Raw signal 50s - spatial average')
% xlabel('time/s')
% ylabel('populational firing rate')
%% Power spectrum
% [N,~]=size(raw_signal_to_process(15000:7000));
% y=fft(raw_signal_to_process(15000:7000));% fft of data
% ps1=abs(y).^2;% power spectrum using fft
% 
% nfft=2^8; %number of fast fourier transforms, the higher this number, the more the frequency resolution of the spectrum goes up
% [ps2,freq2]=pwelch(raw_signal_to_process(15000:7000),chebwin(nfft,100),[],N,Fs);% plotting half of the power spectrum with 50% overlap and chebwin window of length 128
% 
% % freq1=(1:N)*Fs/N;%frequency vector
% figure()
% subplot(3,1,1);
% ps1=ps1(1:length(ps2));
% plot(freq2,20*log(ps1),'b')
% title('POWER SPECTRUM USING FFT METHOD')
% 
% subplot(3,1,2);
% plot(freq2,10*log10(ps2),'r')
% title('POWER SPECTRUM USING PWELCH METHOD-chebyshev window') % Klaus's group usually use this
% [ps3,freq3]=pwelch(raw_signal_to_process,hamming(nfft),[],N,Fs);
% subplot(3,1,3);
% plot(freq3,10*log10(ps3),'r')

nfft            = 4000 ; % for freq resolution = 0.125Hz (fs(1000)/0.125 = 8000) \ Controls freq resolution of psd (most commonly preferred resolution is 0.1Hz)
% raw data to be in Nx1 format; hamming windowed PSD (other methods also available)
[pxx,f]         = pwelch(v1,hamming(nfft),0,nfft,Fs,'psd');

% plot
plot(f,10*log10(pxx),'g','linewidth',1.8,'color',[0.18,0.5,0.1])
xlim([1 40]) % adjust for visualization
legend off
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
% title('Power spectral density') % change according to convienience
%% get amplitude envelope in different bands  input?spike series

% AmplitudeEnvelope_alpha = get_amplitude_envelope(raw_signal_to_process, Fs, alpha_hp, alpha_lp);  % alpha_lp -> hp; alpha_hp -> lp
% AmplitudeEnvelope_beta  = get_amplitude_envelope(raw_signal_to_process, Fs, beta_hp, beta_lp);
% AmplitudeEnvelope_gamma = get_amplitude_envelope(raw_signal_to_process, Fs, gamma_hp, gamma_lp);
AmplitudeEnvelope_delta = get_amplitude_envelope(raw_signal_to_process, Fs, delta_hp, delta_lp);
% AmplitudeEnvelope_theta = get_amplitude_envelope(raw_signal_to_process, Fs, theta_hp, theta_lp);
% 

% % % % alpha_filtered = filter_fir(raw_signal_to_process(:,:),alpha_hp,alpha_lp,Fs,2/alpha_hp);
% figure('Name','signal profile - V1','Position',[50,50,1020,250])
% centered_sig=raw_signal_to_process- mean(raw_signal_to_process);
% norm_cent_sig=raw_signal_to_process/max(abs(raw_signal_to_process));
% centered_env=AmplitudeEnvelope_alpha- mean(AmplitudeEnvelope_alpha);
% AmplitudeEnvelope_alpha=centered_env/max(abs(centered_env));
% plot(time(15000:18000),centered_sig(15000:18000),'linewidth',0.7,'color',[0 0.78 0.55]);
% hold on
% plot(time(15000:18000),AmplitudeEnvelope_alpha(15000:18000),'linewidth',1.2,'color',[0.18 0.5 0.1]);
% % hold on
% 
% % legend(h2,{'alpha','beta','gamma','delta','theta'},'Location','n  orthwest')
% xlabel('time/s')
% ylabel('Average firing rate')
% % title('amplitde envelope - 8Hz-16Hz')
% % 

delta_filtered = filter_fir(raw_signal_to_process(:,:),delta_hp,delta_lp,Fs,2/delta_hp);
figure('Name','signal profile - V1','Position',[50,50,1020,250])
plot(time(15000:28000),delta_filtered(15000:28000),'linewidth',0.7,'color',[0 0.78 0.55]);
hold on
plot(time(15000:28000),AmplitudeEnvelope_delta(15000:28000),'linewidth',1.2,'color',[0.18 0.5 0.1]);
xlabel('time/s')
ylabel('Average firing rate')

% dfa_compute_interval = [dfa_calc_smallest_window dfa_calc_largest_window];
% dfa_fit_interval = [dfa_fit_smallest_window,dfa_fit_largest_window];
% 
% [DFA_exp_alpha,~,~,~] = calculate_DFA(AmplitudeEnvelope_alpha,dfa_compute_interval,dfa_fit_interval,Fs,dfa_calc_window_overlap);
% [DFA_exp_beta,~,~,~]  = calculate_DFA(AmplitudeEnvelope_beta,dfa_compute_interval,dfa_fit_interval,Fs,dfa_calc_window_overlap);
% [DFA_exp_delta,~,~,~] = calculate_DFA(AmplitudeEnvelope_delta,dfa_compute_interval,dfa_fit_interval,Fs,dfa_calc_window_overlap);
% [DFA_exp_theta,~,~,~] = calculate_DFA(AmplitudeEnvelope_theta,dfa_compute_interval,dfa_fit_interval,Fs,dfa_calc_window_overlap);
% [DFA_exp_gamma,~,~,~] = calculate_DFA(AmplitudeEnvelope_gamma,dfa_compute_interval,dfa_fit_interval,Fs,dfa_calc_window_overlap);
% 
% %% Kappa calculation
% if calc_kappa
%     try
%         % avalanche analysis is performed on the raw spike time
%         % series, without any noise added
%         if kappa_threshold_type==1
%             threshold=median(raw_signal_to_process)/2;
%         else
%             threshold=0;
%         end
% 
%         [Kappa_size,Kappa_dur,~] = calculate_kappa(raw_signal_to_process',threshold,kappa_bins,kappa_avalanche_size_min,kappa_avalanche_size_max);
%     catch E
%         Kappa_size=NaN;
%         Kappa_dur=NaN;
%         disp(getReport(E));
%     end
% end