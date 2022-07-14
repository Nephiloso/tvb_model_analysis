function [x_output,y_output]=Logbinning_v2(x,y,res_logbin)
    %% This function transforms the binned data in x,y into log-equidistance
    %% data in x_output and y_output
    %Simon-Shlomo Poil 2006

    %% define the new x axis which should have log-equidistances bin sizes
    d1 = floor(log10(min(x(x>0))));
    d2 = ceil(log10(max(x)));
    x_log_t = logspace(d1,d2(1),(d2-d1)*(res_logbin)+(d2-d1));
    x_log = x_log_t((min(x-(x_log_t(end)-x_log_t(end-1))) <= x_log_t & x_log_t <= max(x+(x_log_t(end)-x_log_t(end-1)))));

    %% Find the midpoint in each bin
    x_log = log10(x_log);
    log_bin_size = x_log(2) - x_log(1);
    midpoint_dist = log_bin_size/2;
    %% define the x_output
    x_output = x_log(1:(end-1))+midpoint_dist;
    x_output = 10.^(x_output);
    %% find the y_output
    x = x(y ~=0);
    y = y(y ~=0);
    m_t = histc(x,10.^(x_log)); %find how many 'x' should be put in each bin
    m = zeros((length(m_t)-1),1);
    m(end) = m_t(end-1) + m_t(end); % remove last point 
    m(1:(end-1)) = m_t(1:(end-2));
    y_index = 0;
    for i=1:(size(m,1))
        if (m(i) == 0)
            y_output(i) = 0; %if no 'x' should be in the bin
        elseif (m(i) == 1)
            y_index = y_index+1; %if 1 'x' should be in the bin
            y_output(i) = y(y_index);
        else
            y_output(i) = 10.^(mean(log10(y((y_index+1):(y_index+m(i)))))); % if more than 1 'x' should be in the bin
        %    y_output(i) = 10.^(median(log10(y((y_index+1):(y_index+m(i)))))); % if more than 1 'x' should be in the bin TEMP!!
            y_index = y_index + m(i);
        end
    end
    return
end