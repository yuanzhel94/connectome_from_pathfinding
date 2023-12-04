function vec2node = vector2node(axon_coord,node_coord)
%vec2node = vector2node(axon_coord,node_coord)
%   compute the vector from each axon to each node
%      Input: 
%        axon_coord: (n_axon,n_dim) matrix of axon location
%        node_coord: (n_node,n_dim) matrix of node location
%      Output:
%        vec2node: (n_node,n_dim,n_axon) matrix of vector (n_dim) from each axon to each node
n_axon = size(axon_coord,1);
n_dim = size(axon_coord,2);
% n_node = size(node_coord,1);
vec2node = node_coord - reshape(axon_coord',1,n_dim,n_axon); % (n_node,n_dim,n_axon)
end