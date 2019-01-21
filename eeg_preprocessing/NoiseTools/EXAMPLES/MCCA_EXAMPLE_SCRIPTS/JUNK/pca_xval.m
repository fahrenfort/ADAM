function [topcs]=nt_pca_xval(x,N,guard);
%[topcs]=nt_pca_xval(x,N,guard) - PCA with cross-validation
%
%  topcs: truncated PCA transform matrix
%  
%  x: data matrix (samplesXchannels or samplesXchannelsXtrials)
%  N: number of chunks
%  guard: number of guard chunks between training data and test chunk

if nargin<3||isempty(guard); guard=1; end
if nargin<2||isempty(N); N=10; end

x=nt_unfold(x);
[nsamples,nchans]=size(x);
nsamples=N*ceil(nsamples/N);
x(nsamples,:)=0; % extends size(x,1) to multiple of N
chunksize=nsamples/N;

% precalculate array of partial covariance matrices
CC=zeros(nchans,nchans,N); % partial covariance matrices
for iChunk=1:N
    xx=x((iChunk-1)*chunksize+(1:chunksize),:);
    CC(:,:,iChunk)=xx'*xx;
end

figure(1); clf

a=[];
for iChunk=1:N
    idx=(iChunk-guard):(iChunk+guard);
    C=sum(CC(:,:,setdiff(1:N,idx)),3);
    [topcs,ev]=nt_pcarot(C);
    aa=diag(topcs'*CC(:,:,iChunk)*topcs);
    aa=flipud(cumsum(flipud(aa)) ./ (1:nchans)');
    a(:,iChunk)=aa/max(aa);
end
C=sum(CC,3); topcs=nt_pcarot(C);
aa=diag(topcs'*C*topcs);
aa=flipud(cumsum(flipud(aa)) ./ (1:nchans)');
aa=aa/max(aa);

semilogy(a); hold on; semilogy(aa, 'linewidth',2);

