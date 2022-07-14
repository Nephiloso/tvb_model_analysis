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

function AggregatePlfData(path)
    files = dir2(path);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    % Print folder names to command window.
    
    crt_folder_path = [path filesep subFolders(1).name];
    simfile_path = [ crt_folder_path filesep 'Simulation.Parameters'];
    
    %check if missing file
    if ~exist(simfile_path,'file') == 2
        ME = MException('Error',...
            'File %s not found',simfile_path);
        throw(ME);
    end
    
    parameters_dict = read_dictionary_from_file( simfile_path );
    
    %these are the parameters that make each simulation unique
    sim_params_string = parameters_dict('Simulation parameters');
    sim_params_array = strsplit(sim_params_string,',');
    
    %get rid of run mode parameter, not useful to distinguish sims
    idx_run_mode = find(strcmp(sim_params_array, 'Run mode'));
    if ~isempty(idx_run_mode)
        sim_params_array = sim_params_array([1:idx_run_mode(1)-1 1:idx_run_mode(1)+1:end]);
    end
    
    num_sims = length(subFolders);
    
    lengthPrePost = 1501;
    
    %amplitude dependence pre-post periods
    table_aDependence = NaN(num_sims,lengthPrePost);
    %plfs for the 10 amplitude percentiels
    table_amplitude_percentile_plfs = NaN(num_sims,10);
    %plfs for bin phases
    table_plfs_phases = NaN(num_sims,33);
    %phase dependence for pre-post periods
    table_pDependence = NaN(num_sims,lengthPrePost);
    %plf for pre-post periods
    table_plfs = NaN(num_sims,lengthPrePost);
    
    table_column_safe_names = strrep(sim_params_array,' ','_');
    
    table_sim_parameters = array2table(NaN(num_sims,length(table_column_safe_names)),'VariableNames',table_column_safe_names);
    
    for k = 1 : num_sims
        crt_folder_path = [path filesep subFolders(k).name];
        simfile_path = [ crt_folder_path filesep 'Simulation.Parameters'];
        
        if ~exist(simfile_path,'file') == 2
            fprintf('Missing simulation parameters file at %s\n',simfile_path);
            continue;
        end
        
        disp(simfile_path);
        
        parameters_dict = read_dictionary_from_file( simfile_path );

        try
            %it has to have external input
            if parameters_dict('External input')==true
                disp(k);
                for p=1:length(sim_params_array)
                    table_sim_parameters{k,p}=parameters_dict(sim_params_array{p});
                end

                load(fullfile(crt_folder_path,'prestim_dependence.mat'),'amplitude_regulation','phase_regulation','plfs','phasebin_plfs','amplitude_percentile_plfs');

                table_aDependence(k,:) = amplitude_regulation;
                table_amplitude_percentile_plfs(k,:) = amplitude_percentile_plfs;
                table_plfs_phases(k,:) = phasebin_plfs;
                table_pDependence(k,:) = phase_regulation;
                table_plfs(k,:) = plfs;
            end
        catch ME
            fprintf('Error while processing directory %s: %s',crt_folder_path,getReport(ME));
        end
        
    end

    save(fullfile(path,'PLF_aggregate.mat'),'-v7.3','table_aDependence','table_amplitude_percentile_plfs','table_plfs_phases','table_pDependence','table_plfs','table_sim_parameters');
end