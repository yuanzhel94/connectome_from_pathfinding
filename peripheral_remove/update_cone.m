function updated_coord = update_cone(cone_coord,unitF,step_length)
%function [update_coord,update_vec] = update_cone(cone_coord,unitF,step_length)
%   Update axons accoding to its previous location, force recieved (unit
%           vector), and step_length.
%   Inputs:
%       cone_coord: (n_axon,n_dim) matrix of axon growth cone coordinates.
%       unitF: (n_dim,n_axon) of the unit force that each axon received
%       step_length: a scalar of step length

%   Outputs:
%       update_coord: (n_axon,n_dim) the new growth cone coordinate
%       (removed) update_vec: the growth vector of each axon, i.e., unitF * step_length

update_vec = unitF * step_length;
updated_coord = cone_coord + update_vec';

end