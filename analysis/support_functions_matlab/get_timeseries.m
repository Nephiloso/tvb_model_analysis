function [raw_network_signal,stimulus_time_series] = get_timeseries(simulation_data_path, simu_file_name, stimu_file_name)
    simulation_data_file_path = fullfile(simulation_data_path,simu_file_name);
    result = readmatrix(simulation_data_file_path);
    v1 = result(:,2);
    raw_network_signal= v1;

    %% Extract stimulus time series
    if contains(simu_file_name,'stim')
        stimulus_time_series = readmatrix(fullfile(simulation_data_path, stimu_file_name));
        stimulus_time_series = floor(stimulus_time_series);

        %turn from 0-based index in python to 1-based index in
        %matlab
        stimulus_time_series = stimulus_time_series + 1;
    else
        %stimulus size is 0, then generate random series of
        %stimuli, and compute the plfs with respect to these
        %times
        no_stims = floor(length(raw_network_signal)/1000)-1;
        stimulusPeriodicity = 1000;
        jitter = 500;
        stimulus_time_series = (1:no_stims)*stimulusPeriodicity - jitter/2 + unidrnd(jitter,1,no_stims);
    end
end