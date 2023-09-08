function [axon_coord,success,axons] = simulate_network_old(r,axon_coord,node_coord,beta,step_length,max_step,self_attract,NameValueArgs)
%[axon_coord,success,axons] = simulate_network(r,axon_coord,node_coord,beta,step_length,max_steps,self_attraction,NameValueArgs)
%   Detailed explanation goes here

% This function has been replaced by simulate_network() for improved computational efficiency.
% Note the "axons" output are different between the two functions

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

nvarg_struct1 = struct; % struct for initial growth
if ~self_attract
    nvarg_struct1.axon_assignment = nearest_node(axon_coord,node_coord);
end

if isfield(NameValueArgs,"weights")
    nvarg_struct1.weights = NameValueArgs.weights;    
end

if isfield(NameValueArgs,"noise")
    noise = NameValueArgs.noise;
    nvarg_struct1.noise_pd = make_noise_pd(noise(1),noise(2));  
end


nvarg_struct2 = nvarg_struct1; % struct for following growth
if isfield(NameValueArgs,"delta_th")
    nvarg_struct2.delta_th = NameValueArgs.delta_th;    
end

nvarg1 = namedargs2cell(nvarg_struct1);

[axon_coord,grow,prev_unitF,axons] = axon_growth_selective_old(r,axons,grow,axon_coord,node_coord,beta,step_length,nvarg1{:});
count = 1;
while any(grow) & (count<max_step)
    nvarg_struct2.prev_unitF = prev_unitF;
    nvarg2 = namedargs2cell(nvarg_struct2);
    [axon_coord,grow,prev_unitF,axons] = axon_growth_selective_old(r,axons,grow,axon_coord,node_coord,beta,step_length,nvarg2{:});
    count = count + 1;
end
success = ~grow;

end