function [c,tw]=nt_xcov(x,y,shifts,w);
%[c,tw]=nt_xcov(x,y,shifts,w) - cross-covariance of X and time-shifted Y
%
%
%  c: cross-covariance matrix
%  tw: total weight
%
%  x,y: data to cross correlate
%  shifts: array of time shifts (must be non-negative)
%  w: weights
%  
% This function calculates, for each pair of columns (Xi,Yj) of X and Y, the
% scalar products between Xi and time-shifted versions of Yj. 
% Shifts are taken from array SHIFTS. 
%
% The weights are applied to X.
%
% X can be 1D, 2D or 3D.  W is 1D (if X is 1D or 2D) or 2D (if X is 3D).
% 
% Output is a 2D matrix with dimensions ncols(X)*(ncols(Y)*nshifts).
%
% NoiseTools

if nargin<4; w=[]; end;
if nargin<3||isempty(shifts); shifts=0; end;

if ~isempty(w) && size(x,1)~=size(w,1); error('X and W should have same nrows'); end
if size(x,3)~=size(y,3); error('X and Y should have same npages'); end
if ~isempty(w) && size(x,3)~=size(w,3); error('X and W should have same npages'); end

shifts=shifts(:); 
nshifts=numel(shifts); 

[mx,nx,ox]=size(x);
[my,ny,oy]=size(y);
c=zeros(nx,ny*nshifts);

if ~isempty(w)
    x=nt_fold(nt_vecmult(nt_unfold(x),nt_unfold(w)),mx);
end

% cross covariance
for k=1:ox
    yy=nt_multishift(y(:,:,k),shifts);
    xx=x(1:size(yy,1),:,k);
    c=c+xx'*yy;
end

if isempty(w)
%     tw=ox*ny*size(yy,1);
    tw=ox*size(yy,1);
else
    w=w(1:size(yy,1),:,:);
    tw=sum(w(:));
end
