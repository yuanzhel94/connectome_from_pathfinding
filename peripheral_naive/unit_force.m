function unitF = unit_force(sumF)
%function unitF = unit_force(sumF)
%   Compute the unit force that each axon received
%   Could also be used to compute the unit vector of axon updates
%       Input:
%           sumF: (n_dim,n_axon) of the sum force that each axon received
%               (from function sum_force())
%       Output:
%           unitF: (n_dim,n_axon) of the unit force that each axon received



magnitude = sqrt(sum(sumF .^ 2,1));

% when no net force received, i.e., magnitude == 0, the column (three dimensions) of axon are all zeros
% choose a random direction to continue (sample vector whose coordinates from standard normal distribution, then normalize the vector to magnitude of 1)
% do not worry about abrupt change in axon direction
% as long as this will be corrected by angle threshold delta_th in function angle_threshold()

zero_ind = ~magnitude;
while any(zero_ind)
    n_dim = size(sumF,1);
    n_zeros = nnz(zero_ind);
    sumF(:,zero_ind) = normrnd(0,1,n_dim,n_zeros);
    magnitude = sqrt(sum(sumF .^ 2,1));
    zero_ind = ~magnitude;
end

unitF = sumF ./ magnitude;


end