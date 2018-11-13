function [y]=nt_xprod2(x1,x2,dsratio)
%[y,ind]=nt_xprod(x1,x2,dsratio) - form all crossproducts
%
%  y: crossproducts 
% 
%  x1,x2: data (time*channels*trials)
%  dsratio: ratio by which to downsample cross-product.

if nargin<3 || isempty(dsratio); dsratio=1; end
if nargin<2 ; error('!'); end

if ndims(x1) ~= ndims(x2); error('!'); end
if size(x1,1)~=size(x2,1); error('!'); end

if ndims(x1)==3
    y=nt_fold(nt_xprod2(nt_unfold(x1),nt_unfold(x2),dsratio),size(x1,1));
else
    [nsamples,nchans1]=size(x1);
    [nsamples,nchans2]=size(x2);
    nsamples=floor(nsamples/dsratio);

    for iChan1=1:nchans1
        xx=bsxfun(@times,x2,x1(:,iChan1));
        y(:,(iChan1-1)*nchans2+(1:nchans2))=nt_dsample(xx,dsratio);
    end
end
