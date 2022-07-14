function dictionary = read_dictionary_from_file( properties_file_path )
    dictionary = containers.Map();
    
    fid = fopen(properties_file_path);
    
    tline = fgets(fid);
    while ischar(tline)
	split_string = strread(strtrim(tline),'%s','delimiter','=');
        if size(strtrim(tline),2)>=0 && tline(1)~='#' && size(split_string,1)>1
            
            key_string = split_string{1};
            value_string = split_string{2};
            
            %try to convert to string
            status=0;
            try
                [converted_num,status]=str2num(value_string);
            catch
            end

            if status==1
                value=converted_num;
            elseif strcmp(value_string,'True')
                value=1;
            elseif strcmp(value_string,'False')
                value=0;
            else 
                value=value_string;
            end
            
            dictionary(strtrim(split_string{1}))=value;
        end
        tline = fgets(fid);
    end
    fclose(fid);
end
