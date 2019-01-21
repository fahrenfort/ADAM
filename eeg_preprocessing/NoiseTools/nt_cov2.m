function [c,tw]=nt_cov2(x,w);
%[c,tw]=nt_cov2(x,w) - weighted covariance
%
%  c: covariance matrix
%  tw: total weight 
%
%  x: data
%  w: weights
%  
% X can be 1D, 2D or 3D.  
% W can be 1D (if X is 1D or 2D) or 2D (if X is 3D). The same weight is
% applied to each column.
% 
% Output is a 2D matrix with dimensions (ncols(X)*numel(SHIFTS))^2.
% It is made up of an ncols(X)*ncols(X) matrix of submatrices, each of 
% dimensions numel(SHIFTS)*numel(SHIFTS).
%
% NoiseTools

if nargin<2; w=[]; end;
if prod(size(x))==0; error('data empty'); end

x=nt_unfold(x);
w=nt_unfold(w);

if isempty(w); w=ones(size(x)); end
if size(w,1)~=size(x,1); error ('!'); end
if size(w,2)==1;
    w=repmat(w,[1,size(x,2)]);
elseif size(w,2)~=size(x,2); 
    error('!');
end

c=zeros(size(x,2));
if isempty(w)
    % no weights
    
    c=x'*x;
    tw=size(x,1)*ones(size(c));
    
else
    % weights
    x=x.*w;
    c=x'*x;
    tw=w'*w;
end

