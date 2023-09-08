function noise_pd = make_noise_pd(noise_std,noise_lim)
%function noise_pd = make_noise_pd(noise_std,noise_lim)

%   This function generate a gaussian noise (mu=0) probability distribution for axon growth
%   Inputs:
%       noise_std: scarlar, the standard deviation of gaussian noise.
%       noise_lim: scaralr, the values that gaussian noise will truncated at, i.e., [-noise_lim,noise_lim].
%   Output:
%       the noise probability distribution noise_pd
noise_pd = truncate(makedist('Normal','mu',0,'sigma',noise_std),-noise_lim,noise_lim);

end