function [squares,quads,D]=nt_qpca(x,npcs,nsmooth,nquads)
%[squares,quads]=nt_qpca(x,npcs,nsmooth,nquads) - quadratic PCA
%
%  squares: linear components closest to largest quadratic component 
%  quads: largest quadratic component(s)
%  D: eigenvalues
%
%  x: data (time*channel*trial)
%  npcs: maximum number of data PCs to use [default: all]
%  nsmooth: square smoothing window to apply to xproducts [default: 1]
%  nquads: number of quadratic components to return [default: 1]
%
% NoiseTools.


if nargin<4||isempty(nquads); nquads=1; end
if nargin<3||isempty(nsmooth); nsmooth=1; end
if nargin<2; error('!'); end
[nsamples,nchans,ntrials]=size(x);

x=[x,ones(nsamples,1,ntrials)*max(abs(x(:)))]; % append a DC component to absorb DC

if nargout==1;
    tosquares=nt_qpca0(x,npcs,nsmooth,nquads);
else
    [tosquares,quads,D]=nt_qpca0(x,npcs,nsmooth,nquads);
    quads=quads(:,2:end,:); % discard first (DC)
end

squares=nt_mmat(x,tosquares);
squares=nt_demean2(squares);




    