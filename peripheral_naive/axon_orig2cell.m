function axons_cell = axon_orig2cell(axons_mat)
%function axons_cell = axon_orig2cell(axons_mat)
%   This function convert a matrix containing axon origin coordinate into a
%       cell array
%   Inputs:
%       axons_mat: (n_axon,n_dim) matrix of axon origin coordinates
%   Output:
%       axons_cell: (n_axon,1) cell array, each cell contains a (1,n_dim)
%           vector of axon origin coordinate

axons_cell = num2cell(axons_mat,2);

end