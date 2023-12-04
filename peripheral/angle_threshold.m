function new_unitF = angle_threshold(unitF,prev_unitF,delta_th)
%function new_unitF = angle_threshold(unitF,prev_unitF,delta_th)
%   This function restrict axon growth angle by a threshold.
%   When the angle between two growth steps <= threshold, nothing happens.
%   When the angle between two growth steps > threshold, make it equal to the threshold.
%   Inputs:
%       unitF: (n_dim,n_axon), for each axon, the unit vector of force describing the present update.
%       prev_unitF: (n_dim,n_axon), for each axon, the unit vector of force describing the previous update
%       delta_th: scalar value in radian of angular threshold between two axon updates.
%
%   Output:
%       new_unitF: (n_dim,n_axon), for each axon, the unit vector of update after thresholding.

dotprod = dot(unitF,prev_unitF);
angle_diff = acos(dotprod); % (1,n_axon) vector of angle difference

flag1 = dotprod <= -1+1e-4; % those in inverse direction (allow for a small difference), not sure turn counter-clockwise or clockwise. Remain unchanged, i.e., follow prev_unitF.
flag2 = (angle_diff > delta_th) & (~flag1);% those alrger than threshold but not equal to pi, turn according to the smaller angle.
% when angle_diff > delta_th, rotate so that new direction follows the
% vector falls in between unitF and prev_unitF, with an angle equals the maximum
% threshold from prev_unitF.

new_unitF = unitF;
% reference: https://stackoverflow.com/questions/22099490/calculate-vector-after-rotating-it-towards-another-by-angle-%CE%B8-in-3d-space
% answer from comingstorm - replace prev_unitF with V, and unitF with D
    % D_perp = D - V * ((D . V)/(V . V))
    % D_perp_scaled = D_perp * (|V|/|D_perp|)
    %result = cos(theta) * V + sin(theta) * D_perp_scaled
D_perp = unitF - prev_unitF .* dotprod;
D_perp_scaled = D_perp ./ sqrt(sum(D_perp.^2,1));
result = cos(delta_th) * prev_unitF + sin(delta_th) * D_perp_scaled;
new_unitF(:,flag2) = result(:,flag2);
new_unitF(:,flag1) = prev_unitF(:,flag1);
end
