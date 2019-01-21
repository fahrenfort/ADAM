function [rms,sd,all]=nt_bsrms(x,N,w)
%[rms,sd,all]=nt_bsrms(x,N,w) - calculate rms, estimate sd using bootstrap
%
%  rms: rms over second dimension of x
%  sd: standard deviation of rms calculated by bootstrap
%  all: matrix of all trials
%  
%  x: matrix of observations (time X repetitions)
%  N: number of bootstrap trials [default: 100]
%  w: weight (time X repetitions)

if nargin <2; N=100; end
if nargin<3; w=[]; end

if ndims(x)>2; 
    x=squeeze(x);
    if ndims(x)>2; 
        error('data must be at most 2D'); 
    end
end
if numel(N)>1; error('!'); end; 

if isempty(w)
    [m,n]=size(x);
    all=zeros(m,N);
    for k=1:N
        idx=ceil(n*rand(1,n));
        all(:,k)=sqrt(mean(x(:,idx).^2,2));
    end
    rms=sqrt(mean(x.^2,2));
    sd=sqrt(mean((all-repmat(rms,1,N)).^2,2));
else
    if size(w,2)==1; w=repmat(w,1,size(x,2)); end
    if size(w)~=size(x); error('!'); end
    [m,n]=size(x);
    all=zeros(m,N);
    for k=1:N
        idx=ceil(n*rand(1,n));
        all(:,k)=sqrt(nt_wmean(x(:,idx).^2,w,(:,idx),2));
    end
    rms=sqrt(nt_wmean(x.^2,w,2));
    sd=sqrt(nt_wmean((all-repmat(rms,1,N)).^2,w,2));
end


