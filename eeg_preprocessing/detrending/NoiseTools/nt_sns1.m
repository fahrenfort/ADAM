function x=nt_sns1(x,nneighbors,skip,w,threshold)
% y=nt_sns1(x,nneigbors,skip,w,threshold) - sensor noise suppression
%
%   y: denoised matrix
%
%   x: matrix  to denoise
%   nneighbors: number of channels to use in projection 
%   skip: number of closest neighbors to skip (default: 0)
%   w : weights (default: all ones)
%   threshold: sharedness threshold (default: 2)
%  
%  This version of SNS first regresses out major shared components.


if nargin<5 || isempty(threshold); threshold=2; end
if nargin<4; w=[]; end
if nargin<3 || isempty(skip); skip=0; end
if nargin<2 || isempty(nneighbors); error('need to specify nneighbors'); end
if ~isempty(w) && sum(w(:))==0; error('weights are all zero!'); end
if ~isempty(find(isnan(x))); error('x contains NANs'); end
if numel(nneighbors)>1 || numel(skip)>1; error('nneighbors and skip must be scalars'); end

xx=nt_pca(nt_normcol(x),[],[],10^-6);   % give each sensor equal weight, PCA
xx=xx(:,find(mean(nt_unfold(xx.^2))>threshold),:); % shared components

if numel(xx)==0; error('!'); end

xxx=nt_tsr(x,xx); % strip data of shared components
clear xx

xxxx=x-xxx;     % shared part
%xxx=nt_sns(nt_sns(nt_sns(xxx,nneighbors,skip,w),nneighbors,skip,w),nneighbors,skip,w); % denoise non-shared part
xxx=nt_sns(xxx,nneighbors,skip,w); % denoise non-shared part
x=xxx+xxxx;       % restore shared part




