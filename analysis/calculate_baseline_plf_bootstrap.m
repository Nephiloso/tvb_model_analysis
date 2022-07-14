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

num_bootstraps = 1000000;
num_stimulations = 1999;

resulting_plfs = zeros(num_bootstraps,1);

trials_in_a_percentile_bin = round(num_stimulations/10);

for b=1:num_bootstraps
    all_trials_phase = rand(trials_in_a_percentile_bin,1)*2*pi-pi;
    resulting_plfs(b) = circ_r(all_trials_phase);
end
sorted_plfs = sort(resulting_plfs);

%95% bin
bin_95 = round(0.95*length(sorted_plfs));

disp(sorted_plfs(bin_95));