function prev_unitF = update_prev_unitF(noisy_coord,cone_coord)
%function prev_unitF = update_prev_unitF(noisy_coord,cone_coord)
%   This function is only needed after a dislocation noise is added
%   Compute the unit vector of this update, serving as prev_unitF in next iteration.
%   Inputs:
%       noisy_coord: (n_axon, n_dim) the location of axon growth cone after adding noise.
%       cone_coord: (n_axon, n_dim) the location of axon growth cone in previous step.
%   Output:
%       prev_unitF: (n_dim, n_axon) unit vector of update direction.

coord_diff = noisy_coord' - cone_coord';
[~,prev_unitF] = unit_force(coord_diff);

end