function [A,y]=nt_whiten(x,N)
%[A,y]=nt_whiten(x,N) - whiten spectrally using pca
%
%  A: whitening matrix (to be applied to time-shifted x)
%  y: whitened signal 
%
%  x: signal to whiten
%  N: order (number of time shifts)
%

if nargin<2; error('!'); end

% calculate covariance across time shifts, looping over columns to save space
sz=size(x);
xx=reshape(x,sz(1),prod(sz(2:end)));
C=zeros(N); % covariance of time-shifted data
for iCol=1:size(xx,2);
    xxx=nt_multishift(xx(:,iCol),0:N-1); 
    C=C+xxx'*xxx;
end
C=C/size(x,1);

% PCA, normalize, inverse
[A,evs]=nt_pcarot(C);
tmp=1./evs; a(find(evs<=0))=0;
A=A*diag(tmp)*pinv(A);

if nargout>1
    % apply whitening matrix, keeping only 1st column
    yy=zeros(size(xx));
    for iCol=1:size(yy,2)
        xxx=nt_multishift(xx(:,iCol),0:N-1);
        yy(1:size(xxx,1),iCol)=xxx*A(:,1);
    end
    y=reshape(yy,[size(yy,1),sz(2:end)]);
end

