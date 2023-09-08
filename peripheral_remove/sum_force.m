function sumF = sum_force(inF)
%function sumF = sum_force(inF)
%   Compute the sumed force that each axon received
%       Input:
%           nodalF: (n_node,n_dim,n_axon) of the received force
%       Output:
%           sumF: (n_dim,n_axon) of the sumed force that each axon received

sumF = reshape(sum(inF,1),size(inF,2),size(inF,3));

end