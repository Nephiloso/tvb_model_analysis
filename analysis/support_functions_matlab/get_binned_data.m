function [x_binned,y_binned] = get_binned_data(x_bins,x,y)
    x_binned = zeros(length(x_bins)-1,1);
    y_binned = zeros(length(x_bins)-1,1);
    
    for bin_idx = 1:length(x_bins)-1
        crt_bin_lower_edge = x_bins(bin_idx);
        crt_bin_upper_edge = x_bins(bin_idx+1);
        x_in_bin = x>=crt_bin_lower_edge & x<crt_bin_upper_edge;
        
        x_binned(bin_idx) = (crt_bin_upper_edge-crt_bin_lower_edge)/2 + crt_bin_lower_edge;
        
        assert(sum(isnan(y))==0,'found nan elements here');
        
        y_binned(bin_idx) = mean(y(x_in_bin));
    end
end

