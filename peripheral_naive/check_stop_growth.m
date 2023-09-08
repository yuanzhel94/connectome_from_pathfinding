function flag = check_stop_growth(cone_coord,r)
%function flag = check_stop_growth(cone_coord,r)
%   This function determines if axons should stop growing based on their
%           updated growth cone coordinate.
%   This function can be used for circle/sphere, but not brain mesh.
%   Inputs: 
%       cone_coord: (n_axon, n_dim) the ccoordinate of growth cones
%       r: the radius of the circle/sphere.
%   Output:
%       flag: (n_axon,1) logical vector, 1 if should stop and 0 otherwise.

% for each cone_coordinate, check their distance d to circle/sphere centre
%   if d>=r, stop growth.
%   allow a small error of 1e-4.

d = sqrt(sum(cone_coord .^ 2,2));
flag = d > r - 1e-4;

end