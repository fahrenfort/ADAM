function [ fitted_data, real_data, sigma, A ] = fit_ctf(CTF)
% function [ fitted_data, real_data, sigma, A ] = fit_ctf(CTF)
% J.J.Fahrenfort, VU, 2016

% what to fit
real_data = CTF;

% set min value to 0 SHOULD I DO THIS?
% real_data = real_data-min(real_data);

% make symmetrical if even number
if mod(numel(real_data),2) == 0
    real_data = [real_data(end) real_data];
end
% mirror
real_data = (real_data(end:-1:1) + real_data)/2;

% normalize SHOULD I DO THIS?
% real_data = real_data/sum(real_data);

% compute
x = 1:numel(real_data);
[sigma,mu,A]=mygaussfit(x,real_data,.3);

% compute fitted real_data
fitted_data = mygausswin(x,sigma,mu,A);

