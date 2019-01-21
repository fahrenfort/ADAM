function [mn,sd,all]=nt_bsmean(x,N,w)
%[mn,sd,all]=nt_bsmean(x,N,w) - calculate mean, estimate sd using bootstrap
%
%  mn: mean of x over second dimension
%  sd: standard deviation from mn of bootstrap trials
%  all: matrix of all bootstrap trials
%  
%  x: matrix of observations (time X repetitions)
%  N: number of bootstrap trials [default: 100]
%  w: weight matrix (time X repetitions)

if nargin <2; N=100; end
if nargin <3; w=[]; end

if ndims(x)>2; 
    x=squeeze(x);
    if ndims(x)>2; 
        error('data must be at most 2D'); 
    end
end
if numel(N)>1; error('!'); end

if isempty(w) % special case for speed
    [m,n]=size(x);
    all=zeros(m,N);
    for k=1:N
        idx=ceil(n*rand(1,n));
        all(:,k)=mean(x(:,idx),2);
    end
    mn=mean(x,2);
    sd=sqrt(mean((all-repmat(mn,1,N)).^2,2));
else
    if size(w,2)==1; w=repmat(w,1,size(x,2)); end
    if size(w) ~= size(x); error('!'); end
    [m,n]=size(x);
    all=zeros(m,N);
    for k=1:N
        idx=ceil(n*rand(1,n));
        all(:,k)=nt_wmean(x(:,idx),w(:,idx),2);
    end
    mn=nt_wmean(x,w,2);
    sd=sqrt( mean((all-repmat(mn,1,N)).^2, 2));
end
    


