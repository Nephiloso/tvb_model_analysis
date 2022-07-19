clear all
% add the analysis folder to path
addpath(genpath('C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\Arthur 2020 code\model code\data_analysis'))

%% Run stimulus response analysis
data_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\optimiz_stim';
figures_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\report\figures';
Stimulus_size = [1; 10; 40; 100; 120; 150; 170];
stim_sizes_name = {'1'; '10'; '40'; '100'; '120'; '150'; '170'};

non_cri_simulation_name = '2022-07-13_c_ee6_c_ei22_stim_size';
cri_simulation_name = '2022-07-05_c_ee11.57199543729056_stim_size';
reg_amp_noncrits=[];
reg_amp_crits=[];
phase_reg_noncrits=[];
phase_reg_crits=[];
color_blue = [0 0.5 1];
color_green = [0.22 1 0.22];
for i=1:length(Stimulus_size)
%plot plf for critical networks at different stimulus sizes
stimulus_size_to_analyze=Stimulus_size(i);
[~,reg_amp_noncrit,~,phase_reg_noncrit] = get_prestim_regulation(data_path,...
    [non_cri_simulation_name,num2str(stimulus_size_to_analyze),'_stim_results.csv'],[non_cri_simulation_name,num2str(stimulus_size_to_analyze),'_tempo.csv']);
[~,reg_amp_crit,~,phase_reg_crit] = get_prestim_regulation(data_path,...
    [cri_simulation_name,num2str(stimulus_size_to_analyze),'_stim_results.csv'],[cri_simulation_name,num2str(stimulus_size_to_analyze),'_tempo.csv']);
reg_amp_noncrits=[reg_amp_noncrits,reg_amp_noncrit];
reg_amp_crits=[reg_amp_crits,reg_amp_crit];
phase_reg_noncrits=[phase_reg_noncrits;phase_reg_noncrit];
phase_reg_crits=[phase_reg_crits;phase_reg_crit];
end

[ps,Ss]=polyfit(Stimulus_size',reg_amp_noncrits,1);
ys = polyval(ps,Stimulus_size);
[pc,Sc]=polyfit(Stimulus_size',reg_amp_crits,1);
yc = polyval(pc,Stimulus_size);
R2s=1 - (Ss.normr/norm(reg_amp_noncrits - mean(reg_amp_noncrits)))^2;
R2c=1 - (Sc.normr/norm(reg_amp_crits - mean(reg_amp_crits)))^2;
sprintf('R2 non-critical:%d',R2s);
sprintf('R2 critical:%d',R2c);
fig=figure('color','w');
plot(Stimulus_size,reg_amp_noncrits,'.','color',color_blue,'MarkerSize',15)
hold on
h1=plot(Stimulus_size,ys,'color',color_blue,'Linewidth',2);
hold on
plot(Stimulus_size,reg_amp_crits,'.','color',color_green,'MarkerSize',15);
hold on
h2=plot(Stimulus_size,yc,'color',color_green,'Linewidth',2);
xlim([0 180])
ylim([-1 1])
set(gca,'fontsize', 16);
xlabel('Stimulation size/neurons');
ylabel('Amplitude regulation');
yticks(linspace(-1,1,11))
legend([h1(1) h2(1)],{'Non-critical','Critical'});
saveas(fig,fullfile(figures_path,'fig_4Fwpx.svg'),'svg');

[ps1,Ss1]=polyfit(Stimulus_size,phase_reg_noncrits,1);
ys1 = polyval(ps1,Stimulus_size);
[pc1,Sc1]=polyfit(Stimulus_size,phase_reg_crits,1);
yc1 = polyval(pc1,Stimulus_size);
R2s1=1 - (Ss1.normr/norm(phase_reg_noncrits - mean(phase_reg_noncrits)))^2;
R2c1=1 - (Sc1.normr/norm(reg_amp_crits - mean(phase_reg_crits)))^2;
sprintf('R2 non-critical:%d',R2s1);
sprintf('R2 critical:%d',R2c1);
fig=figure('color','w');
plot(Stimulus_size,phase_reg_noncrits,'.','color',color_blue,'MarkerSize',15)
hold on
h1=plot(Stimulus_size,ys1,'color',color_blue,'Linewidth',2);
hold on
plot(Stimulus_size,phase_reg_crits,'.','color',color_green,'MarkerSize',15)
hold on
h2=plot(Stimulus_size,yc1,'color',color_green,'Linewidth',2);
xlim([0 180])
ylim([0 0.5])
set(gca,'fontsize', 16);
xlabel('Stimulation size/neurons');
ylabel('Phase regulation');
yticks(linspace(0,0.5,6))
yticklabels({'0','0.1','0.2','0.3','0.4','0.5'});
legend([h1(1) h2(1)],{'Non-critical','Critical'});
saveas(fig,fullfile(figures_path,'fig_5Fwpx.svg'),'svg');

fig=figure('color','w');
x = reallog(Stimulus_size);
h1=plot(x,reg_amp_noncrits,'.-','MarkerSize',15,'color',color_blue,'Linewidth',2.5);
hold on
h2=plot(x,reg_amp_crits,'.-','MarkerSize',15,'color',color_green,'Linewidth',2.5);
ylim([-1 1])
set(gca,'fontsize', 16);
xticks(reallog([1 10 40 180]));
xlabel('Stimulation size/neurons');
ylabel('Amplitude regulation');
yticks(linspace(-1,1,11))
legend([h1(1) h2(1)],{'Non-critical','Critical'});
saveas(fig,fullfile(figures_path,'fig_4Fwpx_log.svg'),'svg');

fig=figure('color','w');
h1=plot(x,phase_reg_noncrits,'.-','MarkerSize',15,'color',color_blue,'Linewidth',2.5);
hold on
h2=plot(x,phase_reg_crits,'.-','MarkerSize',15,'color',color_green,'Linewidth',2.5);
xticks(reallog([1 10 40 180]));
ylim([0 0.5])
set(gca,'fontsize', 16);
xlabel('Stimulation size/neurons');
ylabel('Phase regulation');
yticks(linspace(0,0.5,6))
yticklabels({'0','0.1','0.2','0.3','0.4','0.5'});
legend([h1(1) h2(1)],{'Non-critical','Critical'});
saveas(fig,fullfile(figures_path,'fig_5Fwpx_log.svg'),'svg');