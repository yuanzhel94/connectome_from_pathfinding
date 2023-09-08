function [c,c_diag,axon_indx] = axons2c_old(axons,success,node_coord,directed)
%[c,c_diag,axon_indx] = axons2c(axons,success,node_coord,directed)
%   This function compatible with output from function simulate_network_old()
%           where output "axons" is an (n_axon,1) cell array.
%   Compute connectivity matrix of the neuronal networks simulated by axon propagation
%   Inputs:
%       axons: cell array (n_axon,1) containing the growth trajectory of all axons
%       success: (n_axon,1) logical label of if axons growth succeed, i.e.,
%               ~grow, where grow is output from axon_growth_selective() function.
%       node_coord: (n_node,n_dim) matrix of node coordinates
%       directed: if true, return directed connectivity, else return undirected connectivity.

%   Output:
%       nodes: coordinates of node centroids, same as [x_cen,y_cen]
%       c: diagnal off weighted connecivity matrix, rows are efferent (outgrowth), columns are
%       afferent (received)
%       c_diag: the diagnal of the connectivity matrix
%       axon_indx: N-by-N cell array, each cell is a vector containing all the axon index connecting the two regions 


axon_start = cell2mat(cellfun(@(x) x(1,:),axons,'UniformOutput', false));
axon_end = cell2mat(cellfun(@(x) x(end,:),axons,'UniformOutput', false));

start_node = nearest_node(axon_start,node_coord); %(n_axon,1) of assignment
end_node = nearest_node(axon_end,node_coord);

n_node = size(node_coord,1);
ind = sub2ind([n_node,n_node],start_node,end_node);

axon_indx = cell(n_node);
c = zeros(n_node);
for i = 1:n_node
    for j = 1:n_node
        assign = ind == sub2ind([n_node,n_node],i,j);
        axon_indx{i,j} = find(assign & success);
        c(i,j) = length(axon_indx{i,j});
    end    
end
c_diag = diag(c);
c = c - diag(diag(c));
if ~directed
    c = c + c';
    axon_indx_cp = axon_indx;
    axon_indx_cp(1:(n_node+1):end) = {[]};
    undirected_indx = cellfun( @(x,y) [x;y], axon_indx,axon_indx_cp','UniformOutput',false);
    axon_indx = undirected_indx;
end
end