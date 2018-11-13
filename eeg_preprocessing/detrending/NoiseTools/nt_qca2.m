function [squares,quad,squares2,quad2,D]=nt_qca2(x,npcs,nsmooth)
%[square,quad,square2,quad2,D]=nt_qca(x,npcs,nsmooth) - maximize induced power using quadratic component analysis
%
%  square: linear component closest to most repeatable quadratic component 
%  quad: most reproducible quadratic component
%  square2: linear component closest to least repeatable quadratic component 
%  quad: least reproducible quadratic component
%  D: eigenvalues
%
%  x: data (time*channel*trial)
%  npcs: maximum number of data PCs to use
%  nsmooth: square smoothing window to apply to xproducts [default: 1]
%
%  Usually we are interested in the first component of 'squares'
%  (component with square closest to best quadratic component). 
% 
% NoiseTools.


if nargin<3||isempty(nsmooth); nsmooth=1; end
if nargin<2; error('!'); end
[nsamples,nchans,ntrials]=size(x);

x=[x,ones(nsamples,1,ntrials)*max(abs(x(:)))]; % append a DC component to absorb DC

[tosquares,quad,tosquares2,quad2,D]=nt_qca02(x,npcs,nsmooth);

squares=nt_mmat(x,tosquares);
squares=nt_demean2(squares);
r=nt_repeatability(squares.^2);

%figure(10); plot(r); pause

[dummy,idx]=sort(r,'descend');
squares=squares(:,idx,:);

squares2=nt_mmat(x,tosquares2(:,1));



    