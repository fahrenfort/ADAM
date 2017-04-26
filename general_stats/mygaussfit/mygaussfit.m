function [sigma,mu,A] = mygaussfit(x,y,h)
% [sigma,mu,A] = mygaussfit(x,y,h)
%
% performs a fit to the function
% y = A * exp( -(x-mu)^2 / (2*sigma^2) )
% the fitting occurs by a polyfit on the ln of the data.
%
% h is the threshold which is the fraction from the maximum y height that
% the function uses for fitting. h should be a number between 0-1. if h is
% not supplied, is set to 0 by default, so that all values used for
% fitting are larger than 0.

% threshold
if nargin<3 || isempty(h)
    h=0; % by default takes everything that is larger than 0
end

% cutting
tokeep = y>max(y)*h;
ynew = y(tokeep);
xnew = x(tokeep);

% fitting
ylog=log(ynew);
xlog=xnew;
p=polyfit(xlog,ylog,2);
A2=p(1);
A1=p(2);
A0=p(3);
sigma=sqrt(-1/(2*A2));
mu=A1*sigma^2;
A=exp(A0+mu^2/(2*sigma^2));
sigma = real(sigma);
