function [y,abscissa]=nt_mmx(x,N)
%[y,abscissa]=nt_mmx(x, N) - calculate min-max pairs
%
%  y: array or matrix of min-max pairs
%  
%  x: data 
%  N: target number of pairs (default: 1000 (or size(x,1) if smaller)
%
%  To plot x cheaply: plot(mmx(x)).
%  To get only maxima: y=mmx(x); y(2:2:end);

if nargin<2; N=1000; end

if ndims(x)==3; 
    y=nt_mmx(x(:,1,1)); % to get size
    y=zeros(size(y,1),size(x,2),size(x,3));
    for k=1:size(x,3)
        [y(:,:,k),abscissa]=nt_mmx(x(:,:,k),N);
    end
    return
end

[m,n]=size(x);

if m<2*N
    y=zeros(2*m,n);
    for k=1:n
        xx=x(:,k);
        xx=[xx,xx]';
        xx=xx(:);
        y(:,k)=xx;
    end
else
    N=ceil(m/floor(m/N));
    K=ceil(m/N);
    y=zeros(2*N,n);
    for k=1:n
        xx=x(:,k);
        xx(m:K*N)=xx(m);
        xx=reshape(xx,K,N);
        xx=[min(xx,[],1);max(xx,[],1)];
        xx=xx(:);
        y(:,k)=xx;
    end
end

abscissa=linspace(0, m, 2*N);
    