function [c,directed_axons,undirected_axons] = axons2c(axons,success,node_coord,return_directed,return_undirected)
%function [c,directed_axons,undirected_axons] = axons2c(axons,success,node_coord,return_directed,return_undirected)
%   This function compatible with output from function simulate_network()
%           where output "axons" is an (n_axon,n_max_step) cell array.
%   Compute connectivity matrix of the neuronal networks simulated by axon propagation
%   Inputs:
%       axons: cell array (n_axon,n_max_step) containing the growth trajectory of all axons
%       success: (n_axon,1) logical label of if axons growth succeed (1 for succeed axons)
%       node_coord: (n_node,n_dim) matrix of node coordinates
%       return_directed: if true, return directed median axon trajectory.
%       return_undirected: if true, return undirected median axon trajectory.

%   Output:
%       nodes: coordinates of node centroids, same as [x_cen,y_cen]
%       c: directed weighted connecivity matrix, rows are efferent (outgrowth), columns are
%       afferent (received)
%       directed_axons: (n_node,n_node) full cell array of directed median axon
%               trajectories. Each cell is a (n_step,n_dim) matrix of trajectory.
%       undirected_axons: (n_node,n_node) uppertriangular cell array of undirected median axon trajectories.
%               Each cell is a (n_step,n_dim) matrix of trajectory.

n_axons = size(axons,1);
n_dim = size(node_coord,2);
axon_start = zeros(n_axons,n_dim);
axon_end = zeros(n_axons,n_dim);
axon_length = zeros(n_axons,1);
for i = 1:n_axons
    axoni = vertcat(axons{i,:});
    axon_start(i,:) = axoni(1,:);
    axon_end(i,:) = axoni(end,:);
    if size(axoni,1) == 1
        axon_length(i) = 0;
    else
        axon_length(i) = sum(sqrt(sum((axoni(2:end,:) - axoni(1:end-1,:)).^2,2)));
    end
end


start_node = nearest_node(axon_start,node_coord); %(n_axon,1) of assignment
end_node = nearest_node(axon_end,node_coord);

n_node = size(node_coord,1);
ind = sub2ind([n_node,n_node],start_node,end_node);

c = zeros(n_node);

if return_directed
    directed_axons = cell(n_node);
    if return_undirected
        undirected_axons = cell(n_node);
        for i = 1:n_node
            for j = i:n_node
                if i==j
                    assign = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    c(i,j) = nnz(assign);
                    median_axon = find_median_axon(axons,axon_length,assign);
                    directed_axons{i,j} = median_axon;
                    undirected_axons{i,j} = median_axon;
                else
                    assign1 = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    assign2 = (ind == sub2ind([n_node,n_node],j,i)) & success;
                    assign = assign1 | assign2;
                    c(i,j) = nnz(assign1);
                    c(j,i) = nnz(assign2);

                    median_axon1 = find_median_axon(axons,axon_length,assign1);
                    median_axon2 = find_median_axon(axons,axon_length,assign2);
                    median_axon = find_median_axon(axons,axon_length,assign);
                    
                    directed_axons{i,j} = median_axon1;
                    directed_axons{j,i} = median_axon2;
                    undirected_axons{i,j} = median_axon;
                end
            end
        end
    else
        undirected_axons = nan;
        for i = 1:n_node
            for j = i:n_node
                if i==j
                    assign = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    c(i,j) = nnz(assign);
                    median_axon = find_median_axon(axons,axon_length,assign);
                    directed_axons{i,j} = median_axon;
                else
                    assign1 = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    assign2 = (ind == sub2ind([n_node,n_node],j,i)) & success;
                    c(i,j) = nnz(assign1);
                    c(j,i) = nnz(assign2);

                    median_axon1 = find_median_axon(axons,axon_length,assign1);
                    median_axon2 = find_median_axon(axons,axon_length,assign2);
                    
                    directed_axons{i,j} = median_axon1;
                    directed_axons{j,i} = median_axon2;
                end
            end
        end
    end
else
    directed_axons = nan;
    if return_undirected
        undirected_axons = cell(n_node);
        for i = 1:n_node
            for j = i:n_node
                if i==j
                    assign = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    c(i,j) = nnz(assign);
                    median_axon = find_median_axon(axons,axon_length,assign);
                    undirected_axons{i,j} = median_axon;
                else
                    assign1 = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    assign2 = (ind == sub2ind([n_node,n_node],j,i)) & success;
                    assign = assign1 | assign2;
                    c(i,j) = nnz(assign1);
                    c(j,i) = nnz(assign2);

                    median_axon = find_median_axon(axons,axon_length,assign);
                    undirected_axons{i,j} = median_axon;
                end
            end
        end
    else
        undirected_axons = nan;
        for i = 1:n_node
            for j = i:n_node
                if i==j
                    assign = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    c(i,j) = nnz(assign);
                else
                    assign1 = (ind == sub2ind([n_node,n_node],i,j)) & success;
                    assign2 = (ind == sub2ind([n_node,n_node],j,i)) & success;
                    c(i,j) = nnz(assign1);
                    c(j,i) = nnz(assign2);
                end
            end
        end
    end
end


end