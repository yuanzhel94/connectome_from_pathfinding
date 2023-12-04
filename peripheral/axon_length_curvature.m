function [axon_length,curv] = axon_length_curvature(median_axons,include)
%function [axon_length,curv] = axon_length_curvature(median_axons,include)
%   This function computes the curvature of axons
%   Inputs:
%       median_axons: (n_node,n_node) cell array, each cell contains the trajectory of the median axon as a matrix
%       include: (n_node,n_node) logical matrix, include this connection into analysis if 1
%   Output:
%       axon_length: (n_node,n_node) matrix, the length of each axon
%       curv: (n_node,n_node) matrix, the mean curvature of each axon (0 for axon trajectory <= 2 steps)


axon_length = zeros(size(include));
curv = zeros(size(include));

indx = find(include);
for i = 1:length(indx)
    indxi = indx(i);
    axoni = median_axons{indxi};
    if ~isempty(axoni)
        [L,R,~] = curvature(axoni);
        axon_length(indxi) = L(end);
        if size(axoni,1)>2
            mean_R = mean(R(2:end-1));
            curv(indxi) = 1/mean_R;
        end
    end
end


end