% This is a demo of how we calculating the CC, CPL, SW, and modularity Q
% Normalized by strength/degree preserved networks
% Here we used a network generated from beta = 1; sl = 1 as an example
% to optimize model parameters for individual connectomes:
% repeat for all networks in a landscape, as well as compute the topological measures for empirical connectomes


%% generate model networks for evaluation
addpath(genpath('peripheral_naive')) %peripheral_naive is the model used in the article

rng(3)%for reproducibility

pert = 1;% perturbation parameter rho
beta = 1;% distance decay of attractive force
sl = 1;% increment step length L_s

n_dim = 2; %2-dimensional circle
r = 30;%radius of circle
max_step = 3 .* r ./ sl;%maximum increment steps allowed
n_nodes = 84;%number of nodes (DK atlas in this case)
n_axons = 2e5;%number of axons to simulate
self_attract = true;%allow an axon originated from node i to receive force from node i
delta_th = pi/12;%maximum angular disparity allowed theta

return_directed = false;%not return the directed median trajectory of axons connecting pairs of nodes
return_undirected = false;%not return the undirected median trajectory of axons connecting pairs of nodes

node_rad = sample_circle_rad(n_nodes,pert);%evenly sample nodes followed by perturbation
node_coord = rad2xy(node_rad,r);%convert polar coordinate to cartesian

axon_rad = sample_circle_rad(n_axons,2); %random sample axons, 2nd argument > 1 for random sampling
axon_coord = rad2xy(axon_rad,r);%convert polar coordinate to cartesian

%simulate the axons: axons - trajectory of each axon; success - if the axon succeed
[axons,success] = simulate_network(r,axon_coord,node_coord,beta,sl,max_step,self_attract,"delta_th",delta_th);

%convert axons to connectivity matrix: 
%c - directed connectivity matrix;
[c,~,~] = axons2c(axons,success,node_coord,return_directed,return_undirected);

%makae undirected connectivity matrix
c_und = c + c' - 2 * diag(diag(c));

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

