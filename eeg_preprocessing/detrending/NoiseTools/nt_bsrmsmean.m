function [r,sd,all]=ft_bsrmsmean(x,N)
%[r,sd,all]=ft_bsrmsmean(x,N) - rms over channels of mean over trials, estimate sd using bootstrap
%
%  r: rms of means over trials
%  sd: standard deviation of r calculated by bootstrap
%  all: matrix of all trials
%  
%  x: matrix of observations (time X repetitions or time X 1 X repetitions)
%  N: number of bootstrap trials [default: 100]

if nargin <2; N=100; end
if ndims(x) ~= 3; error('expected data to be 3D'); end

[m,n,o]=size(x);
all=zeros(m,N);
for k=1:N
    idx=ceil(o*rand(1,o));
    all(:,k)=rms(mean(x(:,:,idx),3),2);
end

r=rms(all,2);
sd=sqrt(mean((all-repmat(r,1,N)).^2,2));


