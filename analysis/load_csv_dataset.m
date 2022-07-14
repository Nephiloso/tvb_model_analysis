function [dataset] = load_csv_dataset(csv_path)
    fid = fopen(csv_path);
    var_names = strsplit(fgetl(fid), ';');
    fclose(fid);

    for i=1:length(var_names)
        var_names{i}=strrep(var_names{i},' ','_');
    end

    opts = detectImportOptions(csv_path);
    dataset = readtable(csv_path,opts);
    dataset.Properties.VariableNames = var_names;
end

