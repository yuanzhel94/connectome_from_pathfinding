function median_axon = find_median_axon(axons,axon_length,assign)
%function median_axon = find_median_axon(axons,axon_length,assign)
%   This function find the axon with median steps
%   Inputs:
%       axons: (n_axons,n_max_steps) cell array of axon trajectories
%       axon_length: (n_axons,1) vector of axon_length
%       assign: (n_axons,1) logical vector of if an axon should be included to compute the median
%
%   Outputs:
%       median_axon: (n_step,n_dim) matrix of the median axon trajectory

lengths = axon_length(assign);
[~,median_id] = min(abs(lengths - median(lengths)));
axons_assigned = axons(assign,:);
median_axon = vertcat(axons_assigned{median_id,:});

end