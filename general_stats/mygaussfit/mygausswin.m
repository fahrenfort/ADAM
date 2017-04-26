function y = mygausswin(x,sigmaval,mu,A)
% function y = mygausswin(x,sigma,mu,A)
% computes a gaussian for the points at x, generates 8 points by default
% sigma is the FWHM
% mu is the center of the window
% A is the amplitude at mu
% By default, y is centered around the middle point of x with a sigma of
% 1.5 and an amplitude of 1
% if you want to later normalize the surface of the gaussian to 1, you can
% do so using normy = y/sum(y);
%
% J.J.Fahrenfort, VU 2016

if nargin < 1 || isempty(x)
    x = 1:8;
end
n = numel(x);
if nargin < 4 || isempty(A)
    A = 1;
end
if nargin < 3 || isempty(mu)
    mu = n/2;
end
if nargin < 2 || isempty(sigmaval)
    sigmaval = 1.5;
end
y = A * exp( -(x-mu).^2 / (2*sigmaval.^2) );