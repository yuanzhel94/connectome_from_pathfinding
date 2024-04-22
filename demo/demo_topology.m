% This is a demo of how we calculating the CC, CPL, SW, and modularity Q
% Normalized by strength/degree preserved networks
% Here we used a network generated from beta = 1; sl = 1 as an example
% to optimize model parameters for individual connectomes:
% repeat for all networks in a landscape, as well as compute the topological measures for empirical connectomes


%% generate model networks for evaluation
addpath(genpath('../peripheral')) 

b = 1;% distance decay of attractive force
sl = 1;% increment step length L_s

% generate network with default settings. To change number of brain regions and other settings, 
% see documentation of function connectomee_from_pathfinding()
[c_und,node_coord,directed_axons,undirected_axons] = connectome_from_pathfinding(b,sl,'rng_state',3); % seed for reproduction
% n_nodes = size(node_coord,1);

%% compute topological properties: CC, CPL, SW, modularity Q

c_standard = 2e5; % constant connectivity that will be scaling to
den = 0.1; % network density of evaluation

cc = zeros(2,1);
cpl = zeros(2,1);
Q = zeros(2,1);

% get the network with specified connectivity
c_th = threshold_proportional(c_und,den);
c_count = sum(c_th(:))/2;
ratio = c_standard / c_count;
c_scale = c_th .* ratio;% linearly scale the connectivity
% because of computational cost in generating null networks for large
% number of networks in a landscape, we generate one null network for each network
% variations in both generated network and null network can be captured by
% multiple landscapes (we used 50 in the study).
% If evaluate empirical network for parameter optimization, repeat the analysis 
% for the number of landscapes in generated networks, take the average to
% capture the variation in null networks.
c_null = null_model_und_sign(c_scale);

cc(1) = mean(clustering_coef_wu(c_scale));
cc(2) = mean(clustering_coef_wu(c_null));

d1 = 1 ./ log(c_scale + 1); % take the 1/log(c) as distance matrix for cpl
d2 = 1 ./ log(c_null + 1);
cpl(1) = charpath(distance_wei(d1),0,0);
cpl(2) = charpath(distance_wei(d2),0,0);

[~,Q(1)] = modularity_und(c_scale);
[~,Q(2)] = modularity_und(c_null);

cc_norm = cc(1) / cc(2)
cpl_norm = cpl(1) / cpl(2)
Q_norm = Q(1) / Q(2)
sw = cc_norm / cpl_norm

