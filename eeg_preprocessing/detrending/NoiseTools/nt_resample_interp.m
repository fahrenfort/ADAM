function y=nt_resample_interp(x,R,method);
% y=nt_resample_interp(x,R,method) - resample with arbitrary ratio
%
%  y: resampled data
%
%  x: data to resample (columnwise)
%  R: ratio of new/old sampling rates
%  method: method to give to interp1 [default: 'spline']
%

if nargin<3||isempty(method); method='spline' ; end

Xq=(1 : 1/R : size(x,1))';
y=interp1(1:size(x,1),x,Xq,method);

