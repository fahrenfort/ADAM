function x=nt_demean2(x,w)
%y=nt_demean2(x,w) - remove mean of each row and page
% 
%  w is optional
%
%  if w is a vector with fewer samples than size(x,1), it is interpreted as
%  a vector of indices to be set to 1, the others being set to 0.
%
% NoiseTools


if nargin<2; w=[]; end
if nargin<1; error('!');end

if ~isempty(w) && numel(w)<size(x,1)
    w=w(:);
    % interpret w as array of indices to set to 1
    if min(w)<1 || max(w)>size(x,1); 
        error('w interpreted as indices but values are out of range');
    end
    ww=zeros(size(x,1),1);
    ww(w)=1;
    w=ww;
end

if ndims(x)==4; 
    for k=1:size(x,4);
        if isempty(w);
            x(:,:,:,k)=nt_demean2(x(:,:,:,k));
        else
            if ndims(w)==4; 
                x(:,:,:,k)=nt_demean2(x(:,:,:,k),w(:,:,:,k));
            else
                x(:,:,:,k)=nt_demean2(x(:,:,:,k),w);
            end
        end
    end
    return
end
            
if ~isempty(w)
    if size(w,3)==1 && size(x,3)~=1;
        w=repmat(w,[1,1,size(x,3)]);
    end
    if size(w,3)~=size(x,3)
        error('W should have same npages as X, or else 1');
    end
end

[m,n,o]=size(x);
if isempty(w)
    x=reshape(nt_demean(reshape(x,m,n*o)), [m,n,o]);
else
    w=repmat(w,[1,n,1]);
    x=reshape(nt_demean(reshape(x,m,n*o),reshape(w,m,n*o)),[m,n,o]);
end
    
