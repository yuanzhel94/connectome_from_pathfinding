function rad = sample_circle_rad(n,perturb)
%function rad = sample_circle_rad(n,perturb)
%  Sample points from a perimeter of circle, return coordinate in radian
%   if perturb = 0, even sample
%   if 0<perturb<=1, even sample and perturb
%   if perturb > 1, random select all points (uniformly)

%   Inputs:
%       n: number of points sampled
%       perturb: determines the perturbation, details above
%   Outputs:
%       rad: (n,1) vector of radian values
if (perturb >= 0) & (perturb <= 1) %even sample and perturb
    % place nodes with even distance on the perimeter
    rad = linspace(0,2*pi*(n-1)/n,n)';

    % perturb all nodes uniformly within a distance (maximum half of its area occupied)
    rad_diff = 2*pi/n;
    pert_val = rad_diff * perturb / 2;
    pert = unifrnd(-1,1,size(rad));
    rad = rad + pert * pert_val;
else %random sample
    rad = rand(n,1) * 2 * pi;
end
end