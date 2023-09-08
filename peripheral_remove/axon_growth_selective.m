function [new_axon_coord,new_nvstruct,grow,axons,cannot_grow] = axon_growth_selective(r,axons,grow,axon_coord,node_coord,beta,step_length,NameValueArgs)
%[new_axon_coord,new_grow,new_prev_unitF,new_axons] = axon_growth_selective(r,axons,grow,axon_coord,node_coord,beta,step_length,NameValueArgs)
%   This function receives required information, only grow axons that need to be
%       grown, and update all required information for next growth.
%    Inputs:
%        r: (scalar) value of the circle/sphere radius
%        axons: cell array (n_axon,n_step), each cell contains the trajectory
%           that has been grown, i.e., a (1,n_dim) matrix.
%        grow: (n_axon,1) vector of logical values, 1 denotes axons that grow in this step
%               (not stopped growth in previous steps). When
%               initiate growth, use a logical vector of 1.
%        axon_coord: (n_1_in_grow,n_dim) matrix of axon growth cone coordinate
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
%    Output:
%        new_axon_coord: (n_1_in_new_grow,n_dim) matrix of updated axon growth cone coordinate
%        grow: (n_axon,1) vector of logical values, i.e., updated grow, where 1 denotes axons need
%               to be grown (not stopped growth in previous steps). 
%        new_nvstruct: contains the name value arguments for the next iteration, including updated prev_unitF, axon_asginment, and others.
%        axons: cell array (n_axon,n_step),updated axons.


arguments
    r
    axons
    grow
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

%!!!!!!!!!!!!!!!all removed
% slicing required variables for axons that need to be grown
% args_s = NameValueArgs;
% axon_coord_s = axon_coord(grow,:);
% if isfield(args_s,'axon_assignment')
%     args_s.axon_assignment = NameValueArgs.axon_assignment(grow,:);
% end
% if isfield(args_s,'prev_unitF')
%     args_s.prev_unitF = NameValueArgs.prev_unitF(:,grow);
%     new_prev_unitF = NameValueArgs.prev_unitF;
% end

% initialize the name value argument struct for next iteration (will be returned)
new_nvstruct = NameValueArgs;

% if delta_th present and prev_unitF not present, remove delta_th to avoid
% error message in grow_axons(), i.e., the first growth step do not have prev_unitF
if isfield(NameValueArgs,'delta_th') & (~isfield(NameValueArgs,'prev_unitF'))
    NameValueArgs = rmfield(NameValueArgs,"delta_th");
end

nvargs = namedargs2cell(NameValueArgs); % flatten with {:} as input to axon_growth()


% grow those axons that need to be grown
% updated_coord - (n_1_in_grow,n_dim);
% new_unitF - (n_dim,n1_in_grow)
[updated_coord,new_unitF,cannot_grow] = grow_axons(axon_coord,node_coord,beta,step_length,nvargs{:});
cannot_grow = cannot_grow';


% find the axons that stop after this growth
stop = check_stop_growth(updated_coord,r); %(n_1_in_grow,1) logical, 1 means stopped in this growth

% find the axon ends for those stopped in this growth
p = find_axon_end(axon_coord(stop,:),updated_coord(stop,:),r);%(n_stop_axon,n_dim)


% for those axons that stop, replace update_coord with their intersection to circle/sphere
updated_coord(stop,:) = p;

not_stop = (~stop) & (~cannot_grow) ;

% update the axon growth cone coordinate for next iteration, i.e.,
% new_axon_coord being the axon_coord in next iteration
% new_axon_coord contains the coordinate of axons that does not stop
new_axon_coord = updated_coord(not_stop,:);


% update prev_unitF (only those not stopped) for next iteration
new_nvstruct.prev_unitF = new_unitF(:,not_stop);

% % % % update axon trajectory for nex iteration
% % % axons_grown = axons(grow,:);% find axons grow in this step (previous trajectory)
% % % % cat their trajectory along a third dimension into a matrix (all did not stop previously, thus each is a n_step * n_dim matrix)
% % % M = cat(3,axons_grown{:}); % (n_step,n_dim,n_1_in_grow)
% % % M = cat(1,M,reshape(updated_coord',1,size(M,2),size(M,3)));% add the new trajectory
% % % axons_grown = squeeze(num2cell(M,[1,2]));% convert to cell array of trajectory
% % % new_axons = axons;
% % % new_axons(grow,:) = axons_grown; %update axons for next iteration

% new approach: update axons, where axons is a (n_axon,n_steps) cell array,
% create a new coloumn for those grow at this step

axons(grow,end+1) = num2cell(updated_coord,2);


% update the grow label for next iteration
grow(grow) = not_stop; % updated grow for next iteration



if isfield(new_nvstruct,'axon_assignment')
    new_nvstruct.axon_assignment = NameValueArgs.axon_assignment(not_stop,:);
end


end