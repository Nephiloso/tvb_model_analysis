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

function [ av_size, av_dura ] = get_avalanches( input, threshold )
    ind = find(input>threshold);
    av_size = zeros(1,size(ind,2));
    av_dura = zeros(1,size(ind,2));
    count = 0;
    i=1;
    while i <= size(ind,2)
        start = i;
        i=i+1;
        %find where avalanche ends
        while i <= size(ind,2) && ind(i) == ind(i-1)+1
            i = i+1;
        end
        count = count+1;
        av_size(count) = sum(input(ind(start):ind(i-1)));
        av_dura(count) = i-start;
    end
    av_size = av_size(1:count);
    av_dura = av_dura(1:count);
end