function [axons,success,accident_stop] = simulate_network(r,axon_coord,node_coord,beta,step_length,max_step,self_attract,NameValueArgs)
%[axon_coord,success,axons] = simulate_network(r,axon_coord,node_coord,beta,step_length,max_steps,self_attraction,NameValueArgs)
%   Detailed explanation goes here

%   Outputs:
%       axons: (n_axon, max_step) cell array if axon growth trajectory.
%       success: (n_axon,1) logical matrix labeling if an axon is success.


arguments
    r
    axon_coord
    node_coord
    beta
    step_length
    max_step
    self_attract
    NameValueArgs.weights
    NameValueArgs.delta_th
    NameValueArgs.noise
end

n_axon = size(axon_coord,1);
axons = axon_orig2cell(axon_coord);
grow = ones(n_axon,1,"logical");

nvarg_struct = struct; % struct for initial growth
if ~self_attract
    nvarg_struct.axon_assignment = nearest_node(axon_coord,node_coord);
end

if isfield(NameValueArgs,"weights")
    nvarg_struct.weights = NameValueArgs.weights;    
end

if isfield(NameValueArgs,"noise")
    noise = NameValueArgs.noise;
    nvarg_struct.noise_pd = make_noise_pd(noise(1),noise(2));  
end

if isfield(NameValueArgs,"delta_th")
    nvarg_struct.delta_th = NameValueArgs.delta_th;    
end

% nvarg_struct2 = nvarg_struct1; % struct for following growth




count = 0;
accident_stop = zeros(size(grow),'logical');
while any(grow) & (count<max_step)
    grow_old = grow;
    nvcell = namedargs2cell(nvarg_struct);
    [axon_coord,nvarg_struct,grow,axons,cannot_grow] = axon_growth_selective(r,axons,grow,axon_coord,node_coord,beta,step_length,nvcell{:});
    accident_stop(grow_old) = cannot_grow;
    count = count + 1;
end
success = (~grow) & (~accident_stop);

end