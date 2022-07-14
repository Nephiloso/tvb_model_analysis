function [cluster_beginnings,cluster_endings] = find_clusters(vector)
    %gets the beginning and ending of sections that start and end with 1
    
    cluster_beginnings = find( ( vector(1:end-1)==0 | isnan(vector(1:end-1)) ) & vector(2:end)==1 )+1;
    cluster_endings = find( ( vector(1:end-1)==1 | isnan(vector(1:end-1)) ) & vector(2:end)==0);
    
    if vector(1)==1
        cluster_beginnings = [1 cluster_beginnings];
    end
    
    if vector(end)==1
        cluster_endings = [cluster_endings length(vector)];
    end
end