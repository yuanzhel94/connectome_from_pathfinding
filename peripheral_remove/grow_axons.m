function [updated_coord,new_unitF,cannot_grow] = grow_axons(axon_coord,node_coord,beta,step_length,NameValueArgs)
%function [updated_coord,new_unitF] = grow_axons(axon_coord,node_coord,beta,step_length,NameValueArgs)
%   This function can be used to initiate axons (first growth step) with
%       force gradient, or grow axons (following growth steps).
%   When initiate growth, do not provide variable delta_th and prev_unitF.
%   Optional Name-value pair arguments could be specified to consider:
%       1) axon_assignment: self_attraction avoidance
%       2) weights: weight the attractive force from all nodes by a vector
%       3) delta_th and prev_unitF: specify maximum allowed angle change
%       4) noise_pd: dislocation noise in growth
%    Inputs:
%        axon_coord: (n_axon,n_dim) matrix of axon growth cone coordinate
%        node_coord: (n_node,n_dim) matrix of node coordinate
%        beta: (scalar) powerlaw decay parameter for attractive force
%        step_length: (scalar) value of step length
%        NameValueArgs
%               axon_assignment (optional): (n_axon,1) specify when want to avoid 
%                   attractive force from the nodes that axons originates from.
%               weights (optional): (n_node,1) of the weights of each node (e.g., weighted by nodal area)
%               delta_th (optional): (scalar) threshold of maximum angle change that growth can make
%               prev_unitF(optional): (n_dim,n_axon) of the unit force that each axon
%                   received in the previous step of growth
%               noise_pd(optional): (distribution object) distribution of noise values
%   Output:
%        updated_coord: (n_axon,n_dim) the growth cone coordinate after growth, may
%               outside the geometry
%        new_unitF: (n_dim,n_axon) the unit vector denotes the growth direction of axons
%               in the present step. Serve as prev_unitF in the next growth
arguments
    axon_coord
    node_coord
    beta
    step_length
    NameValueArgs.axon_assignment
    NameValueArgs.weights
    NameValueArgs.delta_th
    NameValueArgs.prev_unitF
    NameValueArgs.noise_pd
end

% compute the force that each axon received from each node
nodalF = force_from_node(axon_coord,node_coord,beta);
if isfield(NameValueArgs,"weights")
    nodalF = weighted_force(nodalF,NameValueArgs.weights);
end
if isfield(NameValueArgs,"axon_assignment")
    nodalF = avoid_self_attraction(nodalF,NameValueArgs.axon_assignment);
end

% compute the unit sum force that each axon received
sumF = sum_force(nodalF);
[zero_force,unitF] = unit_force(sumF);

% optional: threshold angle
dir_invert = zeros(size(zero_force),"logical");
if isfield(NameValueArgs,"delta_th")
    if isfield(NameValueArgs,"prev_unitF")
        [dir_invert,unitF] = angle_threshold(unitF,NameValueArgs.prev_unitF,NameValueArgs.delta_th);
    else
        error("missing argument prev_unitF, i.e., the unit force vectors in the previous step propagation");
    end
end

% update the growth cone coordinate
updated_coord = update_cone(axon_coord,unitF,step_length);

% optional: add a dislocation noise
if isfield(NameValueArgs,"noise_pd")
    updated_coord = add_dislocation_noise(updated_coord,NameValueArgs.noise_pd);
    new_unitF = update_prev_unitF(updated_coord,axon_coord);
else
    new_unitF = unitF;
end

cannot_grow = zero_force | dir_invert;

updated_coord(cannot_grow,:) = axon_coord(cannot_grow,:);

end