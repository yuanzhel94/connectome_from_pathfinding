function f = show_median_axons(c,median_axons,show,n_dim)
%function [f,axon_length] = show_median_axons(c,median_axons,show,n_dim)
%   Detailed explanation goes here
%   Inputs:
%       c: (n_node,n_node) matrix of the weight of each connection.
%       median_axons: (n_node,n_node) cell array of median axon trajectories. Each cell is a (n_step,n_dim) matrix of trajectory.
%       show: (n_node,n_node) logical matrix, show the corresponding median axon if 1 and not show if 0.
%       n_dim: dimension that the trajectory is in.
%   Outputs:
%       f: the figure object of plotted median axons.
%       axon_length: the length of each axon
c = c - diag(diag(c));
max_linewidth = 5;
max_weight = max(max(c));
max_weight_log = log10(max_weight);
log_c = log10(c);

% axon_length = zeros(size(show));

show_indx = find(show);

f = figure;
hold;
title(sprintf('max weight %d',max_weight));

for i = 1:length(show_indx)
    indx = show_indx(i);
    if c(indx) ~= 0
        axoni = median_axons{indx};
%         axon_length(indx) = sum(sqrt(sum((axoni(2:end,:) - axoni(1:end-1,:)).^2,2)));
        weight_log = log_c(indx);
        color_code = [weight_log/max_weight_log,0,0];
        if n_dim == 2
            plot(axoni(:,1),axoni(:,2),'color',color_code,'LineWidth',weight_log/max_weight_log*(max_linewidth-1)+1);
        elseif n_dim == 3
            plot3(axoni(:,1),axoni(:,2),axoni(:,3),'color',color_code,'LineWidth',weight_log/max_weight_log*(max_linewidth-1)+1);
        else
            error('plot in %d dimension not supported',n_dim);
        end
    end
end

axis square;
set(gca,'FontSize',25,'FontWeight','bold');
if n_dim == 3
    view(n_dim);
end
end