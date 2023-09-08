function [c,length_struct,curv_avgR,curv_avgC] = axons2c_more(axons,success,node_coord,min_length,direct)
%function [c,length_struct,curv_avgR,curv_avgC] = axons2c_more(axons,success,node_coord)

%   Compute connectivity matrix of the neuronal networks simulated by axon propagation
%   This is different to axons2c function:
%       This function returns the statistic of all axons connecting node pairs (i.e., mean,median,max,min,std of length and curvature)
%       axons2c returns the axon trajectory with median length

n_axons = size(axons,1);
n_dim = size(node_coord,2);
axon_start = zeros(n_axons,n_dim);
axon_end = zeros(n_axons,n_dim);
axon_length = zeros(n_axons,1);
axon_avgR = zeros(n_axons,1);
axon_avgC = zeros(n_axons,1);
axon_steps = zeros(n_axons,1);

for i = 1:n_axons
    axoni = vertcat(axons{i,:});
    axon_start(i,:) = axoni(1,:);
    axon_end(i,:) = axoni(end,:);
    n_steps = size(axoni,1);
    axon_steps(i) = n_steps;
    if n_steps == 1
        axon_length(i) = 0;

    else
%         axon_length(i) = sum(sqrt(sum((axoni(2:end,:) - axoni(1:end-1,:)).^2,2)));
        [L,R,~] = curvature(axoni);
        axon_length(i) = L(end);
        if n_steps > 2
            R_select = R(2:end-1);
            axon_avgR(i) = 1/mean(R_select);
            axon_avgC(i) = mean(1./R_select);
        end

    end
end

start_node = nearest_node(axon_start,node_coord); %(n_axon,1) of assignment
end_node = nearest_node(axon_end,node_coord);
both_nodes = [start_node,end_node];

n_node = size(node_coord,1);

c = zeros(n_node);

eval_pass = (axon_length>min_length) & (n_steps > 2) & success;

if ~direct
    both_nodes = sort(both_nodes,2);
end
ind = sub2ind([n_node,n_node],both_nodes(:,1),both_nodes(:,2));


mean_curv_avgR = zeros(n_node,n_node);
median_curv_avgR = zeros(n_node,n_node);
max_curv_avgR = zeros(n_node,n_node);
min_curv_avgR = zeros(n_node,n_node);
std_curv_avgR = zeros(n_node,n_node);

mean_curv_avgC = zeros(n_node,n_node);
median_curv_avgC = zeros(n_node,n_node);
max_curv_avgC = zeros(n_node,n_node);
min_curv_avgC = zeros(n_node,n_node);
std_curv_avgC = zeros(n_node,n_node);

mean_length = zeros(n_node,n_node);
median_length = zeros(n_node,n_node);
max_length = zeros(n_node,n_node);
min_length = zeros(n_node,n_node);
std_length = zeros(n_node,n_node);


for i = 1:n_node
    for j = 1:n_node
        if i~=j

            indxk = sub2ind([n_node,n_node],i,j);
            labelk = (ind==indxk) & eval_pass;

            if any(labelk)

                c(i,j) = nnz(labelk);

                axon_avgR_k = axon_avgR(labelk);
                axon_avgC_k = axon_avgC(labelk);

                mean_curv_avgR(i,j) = mean(axon_avgR_k);
                mean_curv_avgC(i,j) = mean(axon_avgC_k);

                median_curv_avgR(i,j) = median(axon_avgR_k);
                median_curv_avgC(i,j) = median(axon_avgC_k);

                max_curv_avgR(i,j) = max(axon_avgR_k);
                max_curv_avgC(i,j) = max(axon_avgC_k);

                min_curv_avgR(i,j) = min(axon_avgR_k);
                min_curv_avgC(i,j) = min(axon_avgC_k);

                std_curv_avgR(i,j) = std(axon_avgR_k);
                std_curv_avgC(i,j) = std(axon_avgC_k);

                lengthk = axon_length(labelk);
                mean_length(i,j) = mean(lengthk);
                median_length(i,j) = median(lengthk);
                max_length(i,j) = max(lengthk);
                min_length(i,j) = min(lengthk);
                std_length(i,j) = std(lengthk);
            end
        end
    end
end

if ~direct
    c = c + c';
    mean_curv_avgR = mean_curv_avgR + mean_curv_avgR';
    median_curv_avgR = median_curv_avgR + median_curv_avgR';
    max_curv_avgR = max_curv_avgR + max_curv_avgR';
    min_curv_avgR = min_curv_avgR + min_curv_avgR';
    std_curv_avgR = std_curv_avgR + std_curv_avgR';

    mean_curv_avgC = mean_curv_avgC + mean_curv_avgC';
    median_curv_avgC = median_curv_avgC + median_curv_avgC';
    max_curv_avgC = max_curv_avgC + max_curv_avgC';
    min_curv_avgC = min_curv_avgC + min_curv_avgC';
    std_curv_avgC = std_curv_avgC + std_curv_avgC';

    mean_length = mean_length + mean_length';
    median_length = median_length + median_length';
    max_length = max_length + max_length';
    min_length = min_length + min_length';
    std_length = std_length + std_length';
end

length_struct = struct;
length_struct.mean = mean_length + mean_length';
length_struct.median = median_length + median_length';
length_struct.max = max_length + max_length';
length_struct.min = min_length + min_length';
length_struct.std = std_length + std_length;

curv_avgR = struct;
curv_avgR.mean = mean_curv_avgR + mean_curv_avgR';
curv_avgR.median = median_curv_avgR + median_curv_avgR';
curv_avgR.max = max_curv_avgR + max_curv_avgR';
curv_avgR.min = min_curv_avgR + min_curv_avgR';
curv_avgR.std = std_curv_avgR + std_curv_avgR';

curv_avgC = struct;
curv_avgC.mean = mean_curv_avgC + mean_curv_avgC';
curv_avgC.median = median_curv_avgC + median_curv_avgC';
curv_avgC.max = max_curv_avgC + max_curv_avgC';
curv_avgC.min = min_curv_avgC + min_curv_avgC';
curv_avgC.std = std_curv_avgC + std_curv_avgC';
end