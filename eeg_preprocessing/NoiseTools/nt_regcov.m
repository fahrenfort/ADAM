function r=nt_regcov(cxy,cyy,keep,threshold)
%r=nt_regcov(cxy,cyy,keep,threshold) - regression matrix from cross covariance
%
%  r: matrix to apply to regressor to best model data
%
%  cxy: cross-covariance matrix between data and regressor
%  cyy: covariance matrix of regressor
%  keep: number of regressor PCs to keep (default: all)
%  threshold: eigenvalue threshold for discarding regressor PCs (default: 0)

if nargin<4; threshold=[]; end
if nargin<3; keep=[]; end
if nargin<2; error('!'); end

% PCA of regressor
[topcs,eigenvalues]=nt_pcarot(cyy);

% discard negligible regressor PCs
if ~isempty(keep)
    keep=max(keep,size(topcs,2));
    topcs=topcs(:,1:keep);
    eigenvalues=eigenvalues(1:keep);
end
if ~isempty(threshold)
    idx=find(eigenvalues/max(eigenvalues)>threshold);
    topcs=topcs(:,idx);
    eigenvalues=eigenvalues(idx);
end

% cross-covariance between data and regressor PCs
cxy=cxy';
r=topcs'*cxy;

% projection matrix from regressor PCs
r=nt_vecmult(r,1./eigenvalues');

% projection matrix from regressors
r=topcs*r;

return




