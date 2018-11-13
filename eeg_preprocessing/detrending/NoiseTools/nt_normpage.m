function y=nt_normpage(x,w)
% y=nt_normpage(x,w) - normalize each page so its weighted msq is 1
% 
%   y: normalized data
%
%   x: data to normalize
%   w: weight [not implemented]
% 
% Weight should be either a column vector, or a matrix (2D or 3D) of same
% size as data.

if nargin<2; w=[]; end

[m,n,o]=size(x);
if isempty(w)
    y=nt_normcol(reshape(x,[m*n,o]));
    y=reshape(y,[m,n,o]);
else
    error('not implemented for weights');
end

    