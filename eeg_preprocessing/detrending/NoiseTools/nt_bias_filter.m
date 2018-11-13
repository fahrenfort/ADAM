function [c0,c1]=nt_bias_filter(x,B,A)
%[c0,c1]=nt_bias_filter(x,B,A) - covariance with and w/o filter bias
%
% x: data 
% B,A: filter coefficients
%
% NoiseTools
% see nt_filter_peak for a second-order resonator

if nargin<3;
    error('!');
end

c0=nt_cov(x);
x=filter(B,A,x);
c1=nt_cov(x);
return



% example:
w=0.3;  % center frequency
Q=10;
bw=f/Q;
[B,A]=secondorderPeak(w,bw);



function [B,A] = secondorderPeak(Wo,BW)
%second order resonator

% Inputs are normalized by pi.
BW = BW*pi;
Wo = Wo*pi;

Ab=(10*log10(2));

Gb   = 10^(-Ab/20);
beta = (Gb/sqrt(1-Gb.^2))*tan(BW/2);
gain = 1/(1+beta);

B  = (1-gain)*[1 0 -1];
A  = [1 -2*gain*cos(Wo) (2*gain-1)];
