function [topcs,pwr,y]=nt_pca0(x,shifts,nkeep,threshold,w)
%[topcs,pwr,y]=nt_pca0(x,shifts,nkeep,threshold,w) - time-shift pca
%
%  topcs: matrix to convert data to PCs
%  pwr: power per PC
%  y: PCs
%
%  x: data matrix
%  shifts: array of shifts to apply
%  nkeep: number of PCs to keep
%  w: weight (see nt_cov)
%  threshold: remove components with normalized eigenvalues smaller than threshold (default: 0)
%
% mean is NOT removed prior to processing


if nargin<1; error('!'); end
if nargin<2||isempty(shifts); shifts=[0]; end
if nargin<3; nkeep=[]; end
if nargin<4||isempty(threshold); threshold=0; end
if nargin<5; w=[]; end

[m,n,o]=size(x);

% remove mean
%x=fold(demean(unfold(x)),size(x,1));

% covariance
if isempty(w);
    c=nt_cov(x,shifts);
else
    c=nt_cov(x,shifts,w);
end

% PCA matrix
if ~isempty(nkeep)
    [topcs,ev]=nt_pcarot(c,nkeep);
else
    [topcs,ev]=nt_pcarot(c);
end

%if ~isempty(nkeep); topcs=topcs(:,1:nkeep); end

% power per PC
pwr=diag(topcs'*c*topcs)/(m*o);
if 0
    idx=find(pwr>=threshold*max(pwr));
    pwr=pwr(idx)';
    topcs=topcs(:,idx);
end

% PCs
if nargout>2
    y=nt_mmat(x,topcs);
end

%% test code
if 0
    x=randn(1000,10);
    [topcs,pwr,y]=nt_pca0(x);
    figure(1); plot(pwr);
    figure(2); subplot 121; plot(y); subplot 122; plot(x*topcs);
end
if 0
    x=zeros(1000,10);
    [topcs,pwr,y]=nt_pca0(x);
    figure(1); plot(pwr);
    figure(2); subplot 121; plot(y); subplot 122; plot(x*topcs);
end
if 0 
    x=sin(2*pi*3*(1:1000)'/1000)*randn(1,10);
    x=2*x + randn(size(x));
     [topcs,pwr,y]=nt_pca0(x);
    figure(1); plot(pwr);
    figure(2); subplot 121; plot(x); subplot 122; plot(x*topcs);
end   



