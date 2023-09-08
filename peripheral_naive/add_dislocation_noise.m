function noisy_coord = add_dislocation_noise(cone_coord,noise_pd)
%function noisy_coord = add_dislocation_noise(cone_coord,noise_std,noise_lim)
%   This function add a gaussian dislocation noise (mu=0) to axon growth
%   Inputs:
%       cone_coord: (n_axon,n_dim) cone coordinate that each axon will grow to
%       noise_pd: the probability distribution of noise from make_noise_pd().

%   Output:
%       noisy_coord: (n_axon,n_dim), the final axon growth cone coordinate after adding the dislocation
%                   noise.

n_axon = size(cone_coord,1);
n_dim = size(cone_coord,2);
noise = random(noise_pd,n_axon,n_dim);
noisy_coord = cone_coord + noise;

end