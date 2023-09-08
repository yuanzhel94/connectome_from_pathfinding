function d2node = distance2node(axon_coord,node_coord)
%d2node = distance2node(axon_coord,node_coord)
%   compute the distance from each axon to each node
%      Input: 
%        axon_coord: (n_axon,n_dim) matrix of axon location
%        node_coord: (n_node,n_dim) matrix of node location
%      Output:
%        d2node: (n_node,1,n_axon) matrix of scalar distance from each axon to each node


%compute the vector from each axon to each node
n_axon = size(axon_coord,1);
n_dim = size(axon_coord,2);
n_node = size(node_coord,1);
vec2node = node_coord - reshape(axon_coord',1,n_dim,n_axon); % ((n_node,n_dim,n_axon)
% identical to vector2node() function

%compute the distance from each axon to each node
d2node = squeeze(sqrt(sum(vec2node .^ 2,2))); % (n_node,n_axon)
d2node = reshape(d2node,n_node,1,n_axon); % (n_node,1,n_axon)

end