function [tosquare,quad,tosquare2,quad2,D]=nt_qca02(x,npcs,nsmooth)
%[tosquare,quad,tosquare2,quad2,D]=nt_qca02(x,npcs,nsmooth) - maximize induced power using quadratic component analysis
%
%  tosquare: matrix to most reproducible induced component
%  quad: most reproducible quadratic component
%  tosquare2: matrix to most reproducible induced component
%  quad2: most reproducible quadratic component
%  D: eigenvalues sorted by absolute value
%
%  x: data (time*channel*trial)
%  npcs: maximum number of data PCs to use [default: all]
%  nsmooth: square smoothing window to apply to xproducts [default: 1]
%  nquads: number of quadratic components to return [default: 1]
% 
% NoiseTools.


if nargin<3; nsmooth=1; end
if nargin<2; npcs=[]; end
[nsamples,nchans,ntrials]=size(x);


% PCA & normalize, select PCs to save space
THRESH=10^-12;
[topcs,pwr]=nt_pca0(x,[],[],THRESH); 
if isempty(npcs) || npcs>size(topcs,2); npcs=size(topcs,2); end
topcs=topcs*diag(1/sqrt(pwr));
topcs=topcs(:,1:npcs); 
x=nt_mmat(x,topcs);

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
c1=xxx'*xxx;

% DSS to find most repeatable
[todss,pwr0,pwr1]=nt_dss0(c0,c1,[],0);
tobest=todss(:,2); % first is DC
toworst=todss(:,end);

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
tosquare=topcs*V;

% form square matrix
ii=1;
A=zeros(npcs);
for k=1:npcs
    for j=1:k
        if j==k;
            A(k,j)=toworst(ii);
        else
            A(k,j)=toworst(ii)/2;
            A(j,k)=toworst(ii)/2;
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
tosquare2=topcs*V;



% best & worst quadratic component(s)
if nargout>1;
    quads=zeros(nsamples-nsmooth+1,1,ntrials);
    quads2=zeros(nsamples-nsmooth+1,1,ntrials);
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
        quad(:,:,k)=nt_mmat(xx,todss(:,1));
        quad2(:,:,k)=nt_mmat(xx,todss(:,end));
    end
end

