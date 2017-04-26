function [ fitted_data, real_data, slope, intercept ] = fit_slope(CTF)
% function [ fitted_data, real_data, slope, intercept ] = fit_slope(CTF)
% J.J.Fahrenfort, VU 2016

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

% take the first half for a linear fit
nEl = round(numel(real_data)/2);
x = 1:nEl;
real_data = real_data(x);
params = polyfit(x,real_data,1);
slope = params(1);
intercept = params(2);

% compute fitted real_data
fitted_data = slope*x + intercept;