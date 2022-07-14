analysis_folder = '/run/media/arthur/SAMSUNG/PhD/projects/CROS versatility simulations/2019.06.14_phasespace_stim/';

all_dirs = dir(analysis_folder);

delete(gcp('nocreate'));
myCluster = parcluster('local');
myCluster.NumWorkers = 6;

parfor i=1:length(all_dirs)
    crt_dir_name = all_dirs(i).name;
    
    disp(crt_dir_name);
    run_analysis(fullfile(analysis_folder,crt_dir_name));
end