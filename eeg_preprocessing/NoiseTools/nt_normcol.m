function [y,norm]=nt_normcol(x,w)
% [y,norm]=nt_normcol(x,w) - normalize each column so its weighted msq is 1
% 
%   y: normalized data
%   norm: vector of norms
%
%   x: data to normalize
%   w: weight
% 
% If x is 3D, pages are concatenated vertically before calculating the
% norm. If x is 4D, apply normcol to each book.
% 
% Weight should be either a column vector, or a matrix (2D or 3D) of same
% size as data.
%
% See nt_normrow, nt_normpage, nt_normpagecol.
%
% NoiseTools

if nargin<2; w=[]; end

if isempty(x); error('empty x'); end

if ndims(x)==4;
    if nargin>1; error('weights not supported for 4D data'); end
    [m,n,o,p]=size(x);
    y=zeros(size(x));
    N=zeros(1,n);
    for k=1:p
    	[y(:,:,:,k),NN]=nt_normcol(x(:,:,:,k));
        N=N+NN.^2;
    end
    return
end

if ndims(x)==3;
    
    % 3D: unfold, apply normcol on 2D, fold
    [m,n,o]=size(x);
    x=nt_unfold(x);
    if isempty(w);
        % no weight 
        [y,NN]=nt_normcol(x);
        N=NN.^2;
        y=nt_fold(y,m);
    else
        % weight
        if size(w,1)~=m; error('weight matrix should have same nrows as data'); end 
        if ndims(w)==2 && size(w,2)==1; 
            w=repmat(w,[1,m,o]);
        end
        if size(w)~=size(w); error('weight should have same size as data'); end
        w=nt_unfold(w);
        [y,NN]=nt_normcol(x,w);
        N=NN.^2;
        y=nt_fold(y,m);
    end

else
    
    % 2D
    [m,n]=size(x);
    if isempty(w)

        % no weight
        %N=sqrt(sum(x.^2)/m);
        %y=vecmult(x,1./N);
        N=(sum(x.^2)/m);
        NN=N.^-0.5;
        NN(find(N==0))=0;
        y=nt_vecmult(x,NN);
       
    else

        % weight
        if size(w,1)~=size(x,1); error('weight matrix should have same ncols as data'); end 
        if ndims(w)==2 && size(w,2)==1; 
            w=repmat(w,1,n);
        end
        if size(w)~=size(w); error('weight should have same size as data'); end
        if size(w,2)==1; w=repmat(w,1,n);end
        %N=sqrt(sum((x.^2).*w)./sum(w));
        %y=vecmult(x,1./N);
        N=(sum((x.^2).*w)./sum(w));
        NN=N.^-0.5;
        NN(find(N==0))=0;
        y=nt_vecmult(x, NN);
        
    end
end

norm=N.^0.5;