function weightedF = weighted_force(nodalF,weights)
%function weightedF = weighted_force(nodalF,weights)
%   Weight the powerlaw decayed attractive force by some weights
%       Inputs:
%           nodalF: (n_node,n_dim,n_axon) of the received force, output from force_from_node()
%           weights: (n_node,1) of the weights of each node (e.g., weighted by nodal area)
%       Outputs:
%           weightedF: nodalF weighted by weights
weightedF = nodalF .* weights;

end