function p = find_axon_end(coord1,coord2,r)
%function p = find_axon_end(coord1,coord2,r)
%   For those axons that stop, find their intersection with circle as axon ends.
%       Inputs:
%           coord1: (n_axon,n_dim) coordinate of the step before axon termination.
%           coord2: (n_axon,n_dim) coordinate of the step of axon termination.
%           r: (scalar) radius of the circle.
%       Output:
%           p: (n_axon,n_dim) coordinate of the intersection/axon ends.

a = sum( (coord2 - coord1) .^ 2, 2);
b = 2 * sum( (coord2 - coord1) .* coord1, 2);
c = sum(coord1 .^ 2,2) - r^2;
u1 = (-b + sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
u2 = (-b - sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);

% u1 and u2 should be: one larger than to 0, one smaller than 0 (P1 not on circle/sphere). 
% positive u correspond to the points of axon ends.
label = u1 > u2; % points that should take u1, otherwise take u2
u = zeros(size(u1));
u(label) = u1(label);
u(~label) = u2(~label);
p = coord1 + u .* (coord2 - coord1);


% Reference: http://paulbourke.net/geometry/circlesphere/#:~:text=Intersection%20of%20a%20Line%20and%20a%20Sphere%20(or%20circle)&text=If%20it%20equals%200%20then,the%20sphere%20at%20two%20points.
% Section: Intersection of a line and a sphere (or circle).
% Let P1(x1,y1,z1) and P2(x2,y2,z2) define the line, P3(x3,y3,z3) define the sphere centre.
%       In our case, P3 = (0,0,0) or (0,0) for sphere and circle.
% Let P =  P1 + u*(P2-P1) define the point of intersection with circle/spehre, where u is a scalar.
% Let:
%       a = (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2;
%       b = 2 * ( (x2-x1)*(x1-x3) + (y2-y1)*(y1-y3) + (z2-z1)*(z1-z3) );
%       c = x3^2 + y3^2 + z3^2 + x1^2 + y1^2 + z1^2 - 2*(x3*x1 + y3*y1 + z3*z1) - r^2;
% The solution becomes:
%       u = ( -b +/- sqrt(b^2 - 4*a*c) ) / (2*a);
% Among the two points, the one that is closer to P2 could be considered as
%       the intersection between our segment P1->P2 and circle/sphere.
end