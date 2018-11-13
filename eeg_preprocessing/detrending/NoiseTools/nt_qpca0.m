function [tosquares,quads,D]=nt_qpca0(x,npcs,nsmooth,nquads)
%[tosquares,quads,D]=nt_qpca0(x,npcs,nsmooth,nquads) - quadratic PCA
%
%  tosquares: matrix to linear components closest to largest quadradric component
%  quads: largest quadratic component(s)
%  D: eigenvalues sorted by absolute value
%
%  x: data (time*channel*trial)
%  npcs: maximum number of data PCs to use [default: all]
%  nsmooth: square smoothing window to apply to xproducts [default: 1]
%  nquads: number of quadratic components to return [default: 1]
% 
% NoiseTools.


if nargin<4; nquads=1; end
if nargin<3; nsmooth=1; end
if nargin<2; npcs=[]; end
[nsamples,nchans,ntrials]=size(x);


% PCA & normalize, select PCs to save space
THRESH=10^-12;
[topcs,pwr]=nt_pca0(x,[],[],THRESH); 
if isempty(npcs) || npcs>size(topcs,2); npcs=size(topcs,2); end
topcs=topcs*diag(1./sqrt(pwr));
topcs=topcs(:,1:npcs); 
x=nt_mmat(x,topcs);

nquads=min(nquads,npcs*(npcs+1)/2-1);

% covariance of cross-products
c0=zeros(npcs*(npcs+1)/2);
xxx=zeros(nsamples-nsmooth+1,npcs*(npcs+1)/2);
for k=1:ntrials
    xx=zeros(nsamples,npcs*(npcs+1)/2);
    ii=1;
    for jj=1:npcs
        for kk=1:jj
            xx(:,ii)=x(:,kk,k).*x(:,jj,k);
            ii=ii+1;
        end
    end
    xx=filter(ones(nsmooth,1)/nsmooth,1,xx);                      % lowpass (demodulate)
    xx=xx(nsmooth:end,:,:);
    xxx=xxx+xx;
    c0=c0+xx'*xx;
end

% DSS to find most repeatable
topcs2=nt_pcarot(c0);
tobest=topcs2(:,2); % first is DC

% form square matrix
ii=1;
A=zeros(npcs);
for k=1:npcs
    for j=1:k
        if j==k;
            A(k,j)=tobest(ii);
        else
            A(k,j)=tobest(ii)/2;
            A(j,k)=tobest(ii)/2;
        end
        ii=ii+1;
    end
end

% eigenvectors & values
[V,D]=eig(A);
D=diag(D);
[dummy,idx]=sort(abs(D),'descend'); 
V=V(:,idx);
D=D(idx);
tosquares=topcs*V;

% answer best quadratic component(s)
if nargout>1;
    quads=zeros(nsamples-nsmooth+1,nquads+1,ntrials);
    for k=1:ntrials
        xx=zeros(nsamples,npcs*(npcs+1)/2);
        ii=1;
        for jj=1:npcs
            for kk=1:jj
                xx(:,ii)=x(:,kk,k).*x(:,jj,k);
                ii=ii+1;
            end
        end
        xx=filter(ones(nsmooth,1)/nsmooth,1,xx);                      % lowpass (demodulate)
        xx=xx(nsmooth:end,:,:);                                       % chop off onset artifact
        quads(:,:,k)=nt_mmat(xx,topcs2(:,1:nquads+1));
    end
end

