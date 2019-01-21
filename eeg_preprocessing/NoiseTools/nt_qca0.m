function [tosquares,quads,D]=nt_qca0(x,npcs,nsmooth,nquads)
%[tosquares,quads,D]=nt_qca0(x,npcs,nsmooth,nquads) - maximize induced power using quadratic component analysis
%
%  tosquares: matrix to most reproducible induced component
%  quads: most reproducible quadratic component(s)
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
topcs=topcs(:,1:npcs); 
x=nt_mmat(x,topcs);

%{
Cross-products are formed trial by trial to save space.
%}

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

% DSS to find most repeatable quadratic form
[todss,pwr0,pwr1]=nt_dss0(c0,c1,[],0);
tobest=todss(:,2); % first is DC

% find linear component with square closest to optimal quad form
[tosquares,D]=nt_quad2square(tobest,'colwise');
tosquares=topcs*tosquares;

% on request, answer best quadratic component(s)
if nargout>1;
    nquads=min(nquads,npcs*(npcs+1)/2-1);
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
        quads(:,:,k)=nt_mmat(xx,todss(:,1:nquads+1));
    end
end

