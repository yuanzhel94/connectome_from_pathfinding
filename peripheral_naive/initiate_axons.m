function cone_coord = initiate_axons(axon_coord,orig_dir,r,step_length)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%   Inputs:
%       axon_coord: (n_axon,n_dim) matrix of axon coordinates. This
%           function is limited to 2D circle/3D sphere axon growth.
if strcmp(orig_dir,'det')
    


elseif strcmp(orig_dir,'rand')
    n_axon = size(axon_coord,1);
    cone_coord = zeros(size(axon_coord));
    good_dir = zeros(n_axon,1,"logical");
    while ~prod(good_dir)
        temp_good_dir = zeros(size(good_dir),"logical"); % initiate to record axons turns to good_dir in this iteration
        indx = find(~good_dir);% find axons still not good_dir
        n_sample = length(indx);% count the number of axons that are still not good_dir
        orig_angle = rand(n_sample,1) * 2 * pi;% for each axon that is still not good_dir, sample a new direction
        diff_coord = [cos(orig_angle), sin(orig_angle)] * step_length;% compute their coordinate difference
        new_coord = axon_coord(indx,:) + diff_coord;% compute their new coordinate after growth
        cone_r = sqrt(sum(new_coord .^ 2,2));% compute their distance to circle/sphere centre
        good_dir(indx) = cone_r < r;% find those locating within circle/sphere, label them as good_dir
        temp_good_dir(indx) = good_dir(indx);% update axons that turns to good_dir in this iteration
        cone_coord(temp_good_dir,:) = new_coord(temp_good_dir(indx),:);% update cone_coord for those labelled as good_dir in this iteration
    end

elseif strcmp(orig_dir,'norm')
    % axon outgrowth norm to its perimeter
    cone_coord = axon_coord - axon_coord ./ sqrt(sum(axon_coord .^ 2,2)) * step_length;
else
    error("invalid orig_dir, choose from 'det','rand','norm'");
end



end

