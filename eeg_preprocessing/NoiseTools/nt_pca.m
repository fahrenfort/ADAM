function [z,idx]=nt_pca(x,shifts,nkeep,threshold,w)
%[z,idx]=nt_pca(x,shifts,nkeep,threshold,w) - time-shift pca
%
%  z: pcs
%  idx: x(idx) maps to z
%
%  x: data matrix
%  shifts: array of shifts to apply
%  keep: number of components shifted regressor PCs to keep (default: all)
%  threshold: discard PCs with eigenvalues below this (default: 0)
%  w: weights
%
% Beware: mean is NOT removed prior to processing.

% TODO: reimplement using nt_pca0
nt_greetings;

if nargin<1; error('!'); end
if nargin<2||isempty(shifts); shifts=[0]; end
if nargin<3; nkeep=[]; end
if nargin<4; threshold=[]; end
if nargin<5; w=[]; end

if isnumeric(x)
    [m,n,o]=size(x);
else
    [m,n]=size(x{1});
    o=length(x);
end


% offset of z relative to x
offset=max(0,-min(shifts));
shifts=shifts+offset;           % adjust shifts to make them nonnegative
idx=offset+(1:m-max(shifts));   % x(idx) maps to z

% % remove mean
% x=nt_fold(nt_demean(nt_unfold(x),w),m);

% covariance
c=nt_cov(x,shifts,w);

% PCA matrix
[topcs,evs]=nt_pcarot(c,nkeep,threshold);

%clf; plot(evs); set(gca,'yscale','log'); pause

% apply PCA matrix to time-shifted data
if isnumeric(x)
    z=zeros(numel(idx),size(topcs,2),o);
    for k=1:o
        z(:,:,k)=nt_multishift(x(:,:,k),shifts)*topcs;
    end
else
    z=[];
    for k=1:o
        z{k}(:,:)=nt_multishift(x{k}(:,:),shifts)*topcs;
    end
end


