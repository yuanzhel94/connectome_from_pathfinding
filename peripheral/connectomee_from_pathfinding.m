function [c_und,node_coord,directed_axons,undirected_axons] = connectomee_from_pathfinding(beta,sl,varargin)
% One-line code that genreates a connectome from pathfinding model
%   [c_und,node_coord,directed_axons,undirected_axons] = connectomee_from_pathfinding(beta,sl,name-value-pairs)
%   Optional name-value pairs include:
%       'r' - radius of circle, default 30
%       'N' - number of nodes, default 84
%       'naxon' - number of axons, default 2e5
%       'rho' - perturbation parameter, default 1
%       'theta' - maximum angle change per step, default pi/12
%       'max_step' - maximum number of growth steps, default 3 * r / sl
%       'return_directed' - return median axon trajectory for directed network, default true
%       'return_undirected' - return median axon trajectory for undirected network, default true
%       'self_attract' - if axons receive attraction from nodes where they origin from, default true
%       'rng_state' - random number generator state, default seed not specified
%   Outputs:
%       c_und - undirected connectivity matrix
%       node_coord - coordinates of nodes
%       directed_axons - median axon trajectory for directed network, NaN if return_directed is false
%       undirected_axons - median axon trajectory for undirected network, NaN if return_undirected is false

pa = inputParser;

% add required arguments
addRequired(pa,'beta');
addRequired(pa,'sl');

% default name-value pairs for independent parameters
default_r = 30; % circle radius
default_N = 84; % # of nodes
default_naxon = 2e5; % # of axons
default_rho = 1; % perturbation parameter rho
default_theta = pi/12; % maximum angle change per step theta
default_self_attract = true;
default_return_directed = true; % return median axon trajectory for directed network
default_return_undirected = true; % return median axon trajectory for undirected network
default_rng = [];

% default name-value pair for dependent parameter - i.e., max growth steps
% initialize as empty and assign value later
default_max_step = [];

% add name-value pairs and parse
addParameter(pa,'r',default_r,@isnumeric);
addParameter(pa,'N',default_N);
addParameter(pa,'naxon',default_naxon);
addParameter(pa,'rho',default_rho);
addParameter(pa,'theta',default_theta);
addParameter(pa,'return_directed',default_return_directed);
addParameter(pa,'return_undirected',default_return_undirected);
addParameter(pa,'max_step',default_max_step);
addParameter(pa,'self_attract',default_self_attract);
addParameter(pa,'rng_state',default_rng);

parse(pa,beta,sl,varargin{:});

% set variables from parser
r = pa.Results.r;
N = pa.Results.N;
naxon = pa.Results.naxon;
rho = pa.Results.rho;
theta = pa.Results.theta;
return_directed = pa.Results.return_directed;
return_undirected = pa.Results.return_undirected;
max_step = pa.Results.max_step;
self_attract = pa.Results.self_attract;

if isempty(pa.Results.max_step)
    max_step = 3 * r / sl;
else
    max_step = pa.Results.max_step;
end

if ~isempty(pa.Results.rng_state)
    rng(pa.Results.rng_state)
end

% generate networks on circle
node_rad = sample_circle_rad(N,rho);%evenly sample nodes followed by perturbation
node_coord = rad2xy(node_rad,r);%convert polar coordinate to cartesian

axon_rad = sample_circle_rad(naxon,2); %random sample axons, 2nd argument > 1 for random sampling
axon_coord = rad2xy(axon_rad,r);%convert polar coordinate to cartesian

%simulate the axons: axons - trajectory of each axon; success - if the axon succeed
[axons,success] = simulate_network(r,axon_coord,node_coord,beta,sl,max_step,self_attract,"delta_th",theta);

%convert axons to connectivity matrix: 
%c - directed connectivity matrix;
%directed_axons: directed median axon trajectory connecting pairs of nodes
%undirected_axons: undirected median axon trajectory connecting pairs of nodes
[c,directed_axons,undirected_axons] = axons2c(axons,success,node_coord,return_directed,return_undirected);

%makae undirected connectivity matrix
c_und = c + c' - 2 * diag(diag(c));

end