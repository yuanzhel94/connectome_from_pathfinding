function outF = avoid_self_attraction(inF,axon_assignment)
%function outF = avoid_self_attraction(inF,axon_assignment)
%   For each axon, remove the attractive force (make it 0) from the node that the axon originates from.
%       Inputs:
%           inF: (n_node,n_dim,n_axon) of the received force, could be output from force_from_node() function (i.e.,nodalF), 
%                       or function weighted_force() (i.e., weightedF).
%           axon_assignment: (n_axon,1) vector of node id, e.g., first axon assigned to node 5, assignment(1) = 5
%                               this could be output from nearest_node() function.
%       Outputs:
%           outF: attractive force where self-attraction values replaced with 0.

[n_node, n_dim, n_axon] = size(inF);
if any(axon_assignment > n_node)
    error("some axons are assigned to nodes that do not exist")
end
sub1 = repmat(axon_assignment',1,n_dim);
sub2 = zeros(1,n_axon*n_dim);
for i = 1:n_dim
    sub2((i-1)*n_axon+1:i*n_axon) = i;
end
sub3 = repmat(1:n_axon,1,n_dim);
ind = sub2ind(size(inF),sub1,sub2,sub3);
outF = inF;
outF(ind) = 0;

end