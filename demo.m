% This is a demo of the model
addpath(genpath('peripheral_naive')) %peripheral_naive is the model used in the article

rng(3)%for reproducibility

tic;
pert = 1;% perturbation parameter rho
beta = 0.98;% distance decay of attractive force
sl = 1;% increment step length L_s

n_dim = 2; %2-dimensional circle
r = 30;%radius of circle
max_step = 3 .* r ./ sl;%maximum increment steps allowed
n_nodes = 84;%number of nodes (DK atlas in this case)
n_axons = 2e5;%number of axons to simulate
self_attract = true;%allow an axon originated from node i to receive force from node i
delta_th = pi/12;%maximum angular disparity allowed theta

return_directed = true;%return the directed median trajectory of axons connecting pairs of nodes
return_undirected = true;%return the undirected median trajectory of axons connecting pairs of nodes

node_rad = sample_circle_rad(n_nodes,pert);%evenly sample nodes followed by perturbation
node_coord = rad2xy(node_rad,r);%convert polar coordinate to cartesian

axon_rad = sample_circle_rad(n_axons,2); %random sample axons, 2nd argument > 1 for random sampling
axon_coord = rad2xy(axon_rad,r);%convert polar coordinate to cartesian

%simulate the axons: axons - trajectory of each axon; success - if the axon succeed
[axons,success] = simulate_network(r,axon_coord,node_coord,beta,sl,max_step,self_attract,"delta_th",delta_th);

%convert axons to connectivity matrix: 
%c - directed connectivity matrix;
%directed_axons: directed median axon trajectory connecting pairs of nodes
%undirected_axons: undirected median axon trajectory connecting pairs of nodes
[c,directed_axons,undirected_axons] = axons2c(axons,success,node_coord,return_directed,return_undirected);

%makae undirected connectivity matrix
c_und = c + c' - 2 * diag(diag(c));

% example - weight-distance associations
fs = 50;%fontsize for figure
d = squareform(pdist(node_coord));
indx = find(c_und~=0);
figure;scatter(d(indx),log10(c_und(indx)),30,[1,1,1] .* 0.5,'filled');
h = lsline;
set(h,'linewidth',5,'color','k');
set(gca,'fontsize',fs,'fontweight','bold');
rval = corr(d(indx),log10(c_und(indx)));
text_dim = [.4,.8,.2,.1];
annotation('textbox',text_dim,'String',sprintf('r = %.2f',rval),'FitBoxToText','on','FontSize',fs,'FontWeight','bold','EdgeColor','none');

%show median axon trajectories
w_max = 1e3;
show_median_axons_specify_max(c_und,undirected_axons,triu(c_und~=0,1),w_max,n_dim);
axis square;axis off;title('');
title('');