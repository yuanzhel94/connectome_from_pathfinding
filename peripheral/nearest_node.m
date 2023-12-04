function nn = nearest_node(axon_coord,node_coord)
%function nn = nearest_node(axon_coord,node_coord)
%   assign each axon coordinate to its nearest node
%      Input: 
%        axon_coord: (n_axon,n_dim) matrix of axon location
%        node_coord: (n_node,n_dim) matrix of node location
%      Output:
%        nn: (n_axon,1) vector of node id, e.g., first axon assigned to node 5, nn(1) = 5

d2node = distance2node(axon_coord,node_coord);%(n_node,1,n_axon) matrix of scalar distance from each axon to each node
[~,nn] = min(d2node,[],1);
nn = squeeze(nn);
end