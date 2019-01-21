function x=nt_normrow(x)
% y=nt_normcol(x) - normalize each row so its msq is 1
% 
%   y: normalized data
%
%   x: data to normalize
%

if ndims(x)>2
    s=size(x);
    x=nt_normrow(reshape(x,[x,prod(s(2:end))]));
    x=reshape(x,s);
else
    w=1./sqrt(mean(x.^2,2));
    w(find(isnan(w)))=0;
    x=bsxfun(@times,x,w);
end

    