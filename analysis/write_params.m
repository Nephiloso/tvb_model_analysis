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

function write_params(myPath,myFname,params_to_write)
    filePath = strcat(myPath,filesep,myFname);
    fid = fopen(filePath,'w');
    try
        all_keys_to_write = keys(params_to_write);
        for i = 1:length(all_keys_to_write)
            crt_key = all_keys_to_write{i};
            crt_value = params_to_write(crt_key);

            if(islogical(crt_value))
                LogicalStr = {'false', 'true'};
                value_string = LogicalStr(crt_value+1);
            elseif isnumeric(crt_value)
                value_string = num2str(crt_value);
            else
                value_string = crt_value;
            end
            fprintf(fid, '%s=%s\n',crt_key,value_string);
        end
    catch E
        disp(getReport(E));
    end

    fclose(fid); 
end
