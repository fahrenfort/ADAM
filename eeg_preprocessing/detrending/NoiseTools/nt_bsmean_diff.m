function [mn,sd]=nt_bsmean_diff(x1,x2,N)
%[mn,sd]=nt_bsmean_diff(x1,x2,N) - calculate mean, estimate sd using bootstrap
%
%  mn: mean of x over second dimension
%  sd: standard deviation from mn of bootstrap trials
%  
%  x1: matrix of observations (time X repetitions)
%  x2: matrix of observations (time X repetitions)
%  N: number of bootstrap trials [default: 100]

if nargin <3; N=100; end
if nargin <2; error('!'); end

if ndims(x1)>2||ndims(x2)>2; error('data must be at most 2D'); end

[m1,n1]=size(x1);
[m2,n2]=size(x2);
if m1~=m2; error('x1 and x2 should have same nrows'); end

all=zeros(m1,N);
for k=1:N
    idx1=ceil(n1*rand(1,n1));
    idx2=ceil(n2*rand(1,n2));
    all(:,k)=mean(x1(:,idx1),2)-mean(x2(:,idx2),2);
end

mn=mean(x1,2)-mean(x2,2);
sd=sqrt(mean((all-repmat(mn,1,N)).^2,2));


