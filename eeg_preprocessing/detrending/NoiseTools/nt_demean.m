function [x,mn]=nt_demean(x,w)
%[y,mn]=nt_demean(x,w) - remove weighted mean over cols
% 
%  w is optional
%
%  if w is a vector with fewer samples than size(x,1), it is interpreted as
%  a vector of indices to be set to 1, the others being set to 0.
%
% NoiseTools

if nargin<2; w=[]; end
if nargin<1; error('!');end
nt_greetings;

if iscell(x);
    cellFlag=1;
    x=nt_cell_to_3D(x); 
    w=nt_cell_to_3D(x);
else
    cellFlag=0;
end

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


if size(w,3)~=size(x,3);
    if size(w,3)==1 && size(x,3)~=1;
        w=repmat(w,[1,1,size(x,3)]);
    else
        error('W should have same npages as X, or else 1');
    end
end

[m,n,o]=size(x);
x=nt_unfold(x);

if isempty(w);
    
    mn=mean(double(x),1);
    x=nt_vecadd(x,-mn);
    
else
    
    w=nt_unfold(w);
    
    if size(w,1)~=size(x,1)
        error('X and W should have same nrows'); 
    end
    
    
    if size(w,2)==1;
        mn=sum(nt_vecmult(double(x),w),1) ./ (sum(w,1)+eps);
    elseif size(w,2)==n;
        mn=sum(x.*w) ./ (sum(w,1)+eps);
    else
        error('W should have same ncols as X, or else 1');
    end

    %y=bsxfun(@minus,x,mn);
    x=nt_vecadd(x,-mn);
    
end

x=nt_fold(x,m);

if cellFlag
    x=nt_3D_to_cell(x);
end
