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

function fig = rasterplot_activity_nostim(data_path,idx_start_short_window)
    simulation_data_file_path = fullfile(data_path,'SimulationData.h5');  
    
    %load spike time series for the excitatory and inhibitory populations
    time_series_exc = h5read(simulation_data_file_path,'/data/main_network_excitatory_spikes_timeseries');
    time_series_inh = h5read(simulation_data_file_path,'/data/main_network_inhibitory_spikes_timeseries');

    %load array of which excitatory neuron spikes
    exc_spike_neurons = h5read(simulation_data_file_path,'/data/main_network_excitatory_spikes_individual_neurons');
    %load array of when the excitatory neuron spikes
    exc_spike_times = h5read(simulation_data_file_path,'/data/main_network_excitatory_spikes_individual_times');
    %load array of which inhibitory neuron spikes
    inh_spike_neurons = h5read(simulation_data_file_path,'/data/main_network_inhibitory_spikes_individual_neurons');
    %load array of when the inhibitory neuron spikes
    inh_spike_times = h5read(simulation_data_file_path,'/data/main_network_inhibitory_spikes_individual_times');
    
    %0-based, turn to 1-based matlab indices
    exc_spike_neurons = exc_spike_neurons + 1;
    inh_spike_neurons = inh_spike_neurons + 1;
    
    %gets which neuron is what type
    types = h5read(simulation_data_file_path,'/networks/main_network/types');
    num_neurons = length(types);
    
    exc_in_types = find(types==1);
    inh_in_types = find(types==0);
    
    exc_spike_neurons = exc_in_types(exc_spike_neurons);
    inh_spike_neurons = inh_in_types(inh_spike_neurons);
    
    %size(in ms) of time window where we display the activity of individual
    %neurons (rasterplot)
    short_window_size = 500;
    short_window_start = idx_start_short_window;
    short_window_end = idx_start_short_window+short_window_size;
    
    %size(in ms) of time interval around the short window where we display
    %total activity
    long_window_size = 5000;
    long_window_start = idx_start_short_window-(long_window_size/2-short_window_size/2);
    long_window_end = idx_start_short_window+(long_window_size/2+short_window_size/2);
    
    %sampling frequency of the signal
    Fs = 1000;

    fig=figure('DefaultAxesFontSize',18);
    
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [2 1 14 6.5]);
    
    %% plot total activity of the excitatory/inhibitory populations over the
    %% long time window
    subplot(4,2,[3 5 7]);
    plot( double(0:long_window_size)/Fs, time_series_exc(long_window_start:long_window_end), 'r');hold on;
    plot( double(0:long_window_size)/Fs, time_series_inh(long_window_start:long_window_end), 'b');
    xlabel('Time (seconds)');
    ylabel('# Spikes');
    xlim([0 long_window_size/Fs]);
    max_firing_rate_inh_exc = max( max(time_series_exc(long_window_start:long_window_end)), max(time_series_inh(long_window_start:long_window_end)));
    ylim([0 max_firing_rate_inh_exc]);
    
    %highlight with transparent rectangle the position of the short time
    %window in the long time window
    rectangle('Position', [(long_window_size/2-short_window_size/2)/Fs, 0, short_window_size/Fs, max_firing_rate_inh_exc], ...
            'Curvature', 0.2, ...
            'FaceColor', [1, 0, 0, 0.1], ...
            'EdgeColor', [1, 0, 0, 0.1]);

    %% plot raster plot of individual neurons activity within the short window
    subplot(4,2,[4 6 8]);
    
    %get spikes in the short time window to make raster plot
    exc_idxes = find(exc_spike_times>=short_window_start & exc_spike_times<=short_window_end);
    inh_idxes = find(inh_spike_times>=short_window_start & inh_spike_times<=short_window_end);
    
    dot_size = 5;
    scatter(double(exc_spike_times(exc_idxes)-idx_start_short_window)/Fs,exc_spike_neurons(exc_idxes),dot_size,'r','filled');hold on;
    scatter(double(inh_spike_times(inh_idxes)-idx_start_short_window)/Fs,inh_spike_neurons(inh_idxes),dot_size,'b','filled');

    xlabel('Time (seconds)');
    ylabel('Neuron');
    
    time_range_seconds = ([0 short_window_size])/Fs;
    xlim(time_range_seconds);
    xticks(time_range_seconds(1):0.1:time_range_seconds(2));
    ylim([1 num_neurons]);
    yticks([0 num_neurons/2 num_neurons]);

    %% plot total activity of excitatory and inhibitory populations within the short window
    subplot(4,2,2);
    idxes_spikes_exc = exc_spike_times(exc_idxes)-idx_start_short_window;
    [unique_time_spikes,~,indexes_time_spikes] = unique(idxes_spikes_exc);
    spike_counts = accumarray(indexes_time_spikes,1);
    spikes_timeseries_short = zeros(501,1);
    spikes_timeseries_short(unique_time_spikes+1)=spike_counts;
    plot(double(0:short_window_size)/Fs,spikes_timeseries_short,'r');hold on;

    idxes_spikes_inh = inh_spike_times(inh_idxes)-idx_start_short_window;
    [unique_time_spikes,~,indexes_time_spikes] = unique(idxes_spikes_inh);
    spike_counts = accumarray(indexes_time_spikes,1);
    spikes_timeseries_short = zeros(501,1);
    spikes_timeseries_short(unique_time_spikes+1)=spike_counts;
    plot(double(0:short_window_size)/Fs,spikes_timeseries_short,'b');hold on;

    yl = ylim;
    ylim([0 max(yl)]);
    xlim([0 short_window_size/Fs]);
    set(gca,'xtick',[]);
    ylabel('# Spikes');
end