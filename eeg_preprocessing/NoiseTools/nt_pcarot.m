function [topcs,eigenvalues]=nt_pcarot(cov,nkeep,threshold,N)
% [topcs,eigenvalues]=pcarot(cov,nkeep,threshold,N) - PCA matrix from covariance
%
%  topcs: PCA rotation matrix
%  eigenvalues: PCA eigenvalues
%  
%  cov: covariance matrix
%  nkeep: number of component to keep
%  thresholds: discard components below this threshold
%  N: eigs' K parameter (if absent: use eig)
%
% NoiseTools

if nargin<4; N=[]; end
if nargin<3; threshold=[]; end
if nargin<2; nkeep=[]; end

if ~isempty(N); 
    [V, S] = eigs(cov,N) ;  
else
    [V, S] = eig(cov) ;  
end

V=real(V);
S=real(S);
[eigenvalues, idx] = sort(diag(S)', 'descend') ;
topcs = V(:,idx);

% truncate
if ~isempty (threshold)
    ii=find(eigenvalues/eigenvalues(1)>threshold);
    topcs=topcs(:,ii);
    eigenvalues=eigenvalues(ii);
end

if ~isempty(nkeep)
    nkeep=min(nkeep,size(topcs,2));
    topcs=topcs(:,1:nkeep);
    eigenvalues=eigenvalues(1:nkeep);
end
