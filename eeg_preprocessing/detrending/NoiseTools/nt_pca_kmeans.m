function [topcs,pwr]=nt_pca_kmeans(x,shifts,nkeep)
%[topcs,pwr]=nt_pca_kmeans(x,nkeep) - PCA preceded by kmeans for speed
%
%  topcs: PCA matrix
%
%  x: data (time*channels or time*channels*trials)
%  nkeep: desired number of PCs
%
% The kmeans implementation of VLFeat appears to be very efficient, 
% so we use it to find a small set of clusters to which we apply PCA.
%
% This differs from the usual way of combining PCA and kmeans, which is
% to use PCA to initialize kmeans.

if nargin<3; error('!'); end
if isempty(shifts); shifts=0; end
if shifts
    x=nt_multishift(x,shifts);
end
if nkeep>size(x,2); error('!'); end

x=nt_unfold(x);
[nsamples,nchans]=size(x);

if nkeep>nsamples; error('nkeep greater than number of samples'); end
if nkeep>nchans; error('nkeep greater than number of channels'); end


% normalize channels
nrm=sqrt(mean(x.^2));
tonrm=diag(1./nrm); % normalization matrix
tonrm(find(isnan(x)))=0;
x=x*tonrm;

% cluster
nclusters=round(sqrt(nchans)); 
[C,A]=vl_kmeans(x,nclusters,'algorithm','elkan');

% perform PCA on each cluster
topcs1=zeros(nchans);
pwr=zeros(1,nchans);
idx=0;
for k=1:nclusters
    %[nclusters k numel(find(A==k))]
    [m,p]=nt_pca0(x(:,find(A==k)));
    topcs1(find(A==k),idx+(1:size(m,2))) = m;
    pwr(idx+(1:size(m,2)))=p;
    idx=idx+size(m,2);
end

% sort by decreasing power
[~,idx]=sort(pwr,'descend');

topcs2=nt_pca0(nt_mmat(x,topcs1(:,idx(1:nkeep))));

topcs=tonrm*(topcs1(:,idx(1:nkeep))*topcs2);








