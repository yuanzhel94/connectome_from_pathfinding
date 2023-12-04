function [axons,success] = simulate_random_walk_network(r,axon_coord,sl,max_step,max_angle)
%[axons,success] = simulate_random_walk_network(r,axon_coord,node_coord,beta,step_length,max_step,max_angle)
%   This function simulate a network where axons are guided by random walk.
%   Max angle is specified to constrain angle change

n_axons = size(axon_coord,1);

success = zeros(n_axons,1,"logical");
axons = cell(n_axons,max_step);

for i = 1:n_axons
    n_step = 1;
    axoni = axon_coord(i,:);
    axons{i,n_step} = axoni;
%     start_node(i) = nearest_node(axoni,node_coord);
    step_rad = rand * 2 * pi;
    new_xy = axoni + rad2xy(step_rad,sl);
    while check_stop_growth(new_xy,r)
        step_rad = rand * 2 * pi;
        new_xy = axoni + rad2xy(step_rad,sl);
    end
    n_step = n_step + 1;
    axoni = cat(1,axoni,new_xy);
    axons{i,n_step} = new_xy;
    while (n_step <= max_step)
        old_rad = step_rad;
        new_rad = unifrnd(-max_angle,max_angle);
        step_rad = old_rad + new_rad;
        new_xy = axoni(end,:) + rad2xy(step_rad,sl);
        if check_stop_growth(new_xy,r)
            new_xy = find_axon_end(axoni(end,:),new_xy,r);
            axoni = cat(1,axoni,new_xy);
            n_step = n_step + 1;
            success(i) = true;
            axons{i,n_step} = new_xy;
            break;
        else
            axoni = cat(1,axoni,new_xy);
            n_step = n_step + 1;
            axons{i,n_step} = new_xy;
        end
    end
%     axon_steps(i) = n_step;
%     if size(axoni,1) > 3
%         axon_avail(i) = true;
%         [L,R,~] = curvature(axoni);
%         axon_length(i) = L(end);
%         R_mean = mean(R(2:end-1));
%         axon_curv_avgR(i) = 1/R_mean;
%         axon_curv_avgC(i) = mean(1./R(2:end-1));
%     end
end



end