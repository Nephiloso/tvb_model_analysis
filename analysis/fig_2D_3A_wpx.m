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

addpath('support_functions_matlab');

data_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\data\optimiz_stim';
figures_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\report\figures';
model_analysis_code_path = 'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes\Arthur 2020 code\model code\data_analysis';

addpath(genpath(model_analysis_code_path));

noncrit_run_number	= 0;

crit_run_number = 1;

stimulus_size_to_analyze = 10;

non_cri_simulation_name = '2022-07-13_c_ee6_c_ei22_stim_size';
cri_simulation_name = '2022-07-05_c_ee11.57199543729056_stim_size';
%plot plf for critical networks at different stimulus sizes
[amplitude_bins_plf_noncrit,reg_amp_noncrit,phase_bins_plf_noncrit,phase_reg_noncrit] = get_prestim_regulation(data_path,...
    [non_cri_simulation_name,num2str(stimulus_size_to_analyze),'_stim_results.csv'],[non_cri_simulation_name,num2str(stimulus_size_to_analyze),'_tempo.csv']);
% [amplitude_bins_plf_crit,reg_amp_crit,phase_bins_plf_crit,phase_reg_crit] = get_prestim_regulation(fullfile(data_path,'runs',...
%     'without individual spikes',sprintf('EC%.2f.IC%.2f.Stim%i.Run%i',crit_exc_conn,crit_inh_conn,stimulus_size_to_analyze,crit_run_number)));
[amplitude_bins_plf_crit,reg_amp_crit,phase_bins_plf_crit,phase_reg_crit] = get_prestim_regulation(data_path,...
    [cri_simulation_name,num2str(stimulus_size_to_analyze),'_stim_results.csv'],[cri_simulation_name,num2str(stimulus_size_to_analyze),'_tempo.csv']);

color_blue = [0 0.5 1];
color_green = [0.22 1 0.22];
% color_red = [1 0.22 0.22];

%plot prestimulus amplitude regulation
% fig = figure('color','w');
% set(fig,'Position',[0 0 600 500]);
% x_percentiles = 10:10:100;
% plot(x_percentiles,amplitude_bins_plf_noncrit,'.','MarkerSize',15,'Color',color_blue);
% h1=lsline;set(h1(1),'color',color_blue,'LineWidth',2);
% hold on;
% plot(x_percentiles,amplitude_bins_plf_crit,'.','MarkerSize',15,'Color',color_green);
% h2=lsline;set(h2(1),'color',color_green,'LineWidth',2);
% % plot(x_percentiles,amplitude_bins_plf_crit,'.','MarkerSize',15,'Color',color_red);
% % h3=lsline;set(h3(1),'color',color_red,'LineWidth',2);
% legend([h1(1) h2(1)],{num2str(round(reg_amp_noncrit,2)),num2str(round(reg_amp_crit,2)),num2str(round(reg_amp_crit,2))});
% xticks(10:10:100);
% ylim([0 0.5]);
% yticks([0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5]);
% yticklabels({'0','','0.1','','0.2','','0.3','','0.4','','0.5'});
% xticklabels({'','','','','50','','','','','100'});
% set(gca,'fontsize', 16);
% xlabel({'Pre-stimulus Amplitude';'(percentile)'});
% ylabel('PLF at 150ms');
% saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_2D.fig']),'fig');
% saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_2D.svg']),'svg');
% [ps,Ss] = polyfit(x_percentiles',amplitude_bins_plf_noncrit,1);
% [pc,Sc] = polyfit(x_percentiles',amplitude_bins_plf_crit,1);

%plot prestimulus phase regulation
fig = figure('color','w');
set(fig,'Position',[0 0 600 500]);
set(gca,'fontsize', 16);
x_phases = -pi:pi/16:pi;
plot(x_phases,phase_bins_plf_noncrit,'.-','MarkerSize',15,'Color',color_blue,'Linewidth',2);
hold on;
plot(x_phases,phase_bins_plf_crit,'.-','MarkerSize',15,'Color',color_green,'Linewidth',2);
% plot(x_phases,phase_bins_plf_crit,'.-','MarkerSize',15,'Color',color_red,'Linewidth',2);
legend({['Noncritical: ',num2str(round(phase_reg_noncrit,2))],['Critical: ',num2str(round(phase_reg_crit,2))]});
xticks(-pi:pi/16:pi);
ylim([0 1]);
yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
yticklabels({'0','','0.2','','0.4','','0.6','','0.8','','1.0'});
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','','0','','\pi'});
xlabel('Phase angle at -5 ms');
ylabel('PLF');
saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_3A.fig']),'fig');
saveas(fig,fullfile(figures_path,['fig_',num2str(stimulus_size_to_analyze),'_3A.svg']),'svg');
