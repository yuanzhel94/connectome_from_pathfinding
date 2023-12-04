function nodalF = force_from_node(axon_coord,node_coord,beta)
%function nodalF = force_from_node(axon_coord,node_coord,beta)
%   Compute the force that each axon received from each node
%      Input: 
%        axon_coord: (n_axon,n_dim) matrix of axon location
%        node_coord: (n_node,n_dim) matrix of node location
%        beta: powerlaw decay parameter for attractive force
%      Output:
%        nodalF: (n_node,n_dim,n_axon) of the received force

%compute the vector from each axon to each node
n_axon = size(axon_coord,1);
n_dim = size(axon_coord,2);
n_node = size(node_coord,1);
vec2node = node_coord - reshape(axon_coord',1,n_dim,n_axon); % ((n_node,n_dim,n_axon)
% identical to vector2node() function

%compute the distance from each axon to each node
d2node = squeeze(sqrt(sum(vec2node .^ 2,2))); % (n_node,n_axon)
d2node = reshape(d2node,n_node,1,n_axon); % (n_node,1,n_axon)
% identical to distance2node() function

%compute the force that each axon reeived from each node
nodalF = vec2node ./ (d2node .^ (beta + 1)); %(n_node,n_dim,n_axon)

end