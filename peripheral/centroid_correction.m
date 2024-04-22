function corrected_rad = centroid_correction(node_rad)
%function corrected_rad = centroid_correction(node_rad)
%   once node coordinates perturbed, the sampled coordinates are not node centroids anymore
%   optionally correct the coordinates

d_clock = node_rad - [node_rad(end)-2*pi;node_rad(1:end-1)];
d_anticlock = [node_rad(2:end); node_rad(1)+2*pi] - node_rad;
adjust = (d_anticlock - d_clock)/4;
corrected_rad = node_rad + adjust;
end