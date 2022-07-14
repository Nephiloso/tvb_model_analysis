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

addpath('CircStat');

num_bootstraps = 10000;
num_stimulations = 1999;

bins = 33;
binsPhaseDep = -pi:pi/16:pi;

resulting_plfs = zeros(num_bootstraps,1);

crt_plfs = zeros(10,1);

binWidth = pi/4;

for b=1:num_bootstraps
    disp(b);
    for t=1:10
        plfs_over_phases = rand(33,1);
        %so here what I say is that the value of phase post stimulus is
        %completely random across trials
        %these makes values for all trials come in the range of -pi/pi
        all_trials_prestim = rand(num_stimulations,1)*2*pi-pi;
        all_trials_poststim = rand(num_stimulations,1)*2*pi-pi;
        plfs_over_phases_binned = zeros(33,1);
        
        for i = 1:length(binsPhaseDep)
            orig = abs(circ_dist(all_trials_prestim,binsPhaseDep(i)));
            pcks = orig < binWidth/2;
            noPcks(i) = nnz(pcks);
        end
        minPcks = min(noPcks);
        
        for i = 1:length(binsPhaseDep)
            orig = abs(circ_dist(all_trials_prestim,binsPhaseDep(i)));
            pcks = orig < binWidth/2;
            pcks = find(pcks);
            pcks = pcks(1:minPcks);
            plfs_over_phases_binned(i)=circ_r(all_trials_poststim(pcks));
        end
        
        crt_plfs(t) = circ_r(binsPhaseDep',plfs_over_phases_binned);
    end
    
    resulting_plfs(b)=mean(crt_plfs);
end
sorted_plfs = sort(resulting_plfs);

%95% bin
bin_95 = round(0.95*length(sorted_plfs));

disp(sorted_plfs(bin_95));