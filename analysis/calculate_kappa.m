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

function [kappa_size,kappa_dur,num_avalanches] = calculate_kappa(data, threshold, noBins, minX, maxX)
    [aval_size, aval_dur] = get_avalanches(data,threshold);
    
    %compute kappa for the avalanche size distribution
    avalanche_size_min = minX;
    avalanche_size_max = maxX;
    avs_with_size_threshold = aval_size(aval_size >= avalanche_size_min & aval_size <= avalanche_size_max);
    t = 10.^(linspace(log10(avalanche_size_min),log10(avalanche_size_max),noBins));
	theoreticalDistribution = zeros(noBins,1);
    actualDistribution = zeros(noBins,1);
    for i = 1:noBins
        theoreticalDistribution(i) = (1 - (avalanche_size_min/avalanche_size_max)^0.5)^-1 * (1 - (minX/t(i))^0.5);
    end
    for i = 1:noBins
        actualDistribution(i) = nnz(avs_with_size_threshold<t(i));
    end
    actualDistribution = actualDistribution / length(avs_with_size_threshold);

    kappa_size = 1 + sum(theoreticalDistribution - actualDistribution)/noBins;
    
    %compute kappa for avalanche duration distribution
    avalanche_duration_min = 1;
    avalanche_duration_max = 500;

    avd_with_duration_threshold = aval_dur(aval_dur >= avalanche_duration_min & aval_dur <= avalanche_duration_max);
    t = 10.^(linspace(log10(minX),log10(maxX),noBins));
    theoreticalDistribution = zeros(noBins,1);
    actualDistribution = zeros(noBins,1);
    for i = 1:noBins
        theoreticalDistribution(i) = (1 - (minX/maxX)^1)^-1 * (1 - (minX/t(i))^1);
    end
    for i = 1:noBins
        actualDistribution(i) = nnz(avd_with_duration_threshold<t(i));
    end
    actualDistribution = actualDistribution / length(avd_with_duration_threshold);
    kappa_dur = 1+ sum(theoreticalDistribution - actualDistribution)/noBins;
    
    num_avalanches = length(aval_size);
end