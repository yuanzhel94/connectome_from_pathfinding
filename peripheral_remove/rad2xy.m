function xy = rad2xy(rad,r)
%function xy = rad2xy(rad,r)
%   convert radian to [x,y] value on a circle of radius r
%   Inputs:
%       rad: (n,1) vector of radian values
%       r: radius of the circle
%   Ouput:
%       xy: (n,2) matrix, x (col1) and y (col2) coordinates of sampled points
xy = [cos(rad)*r,sin(rad)*r];
end