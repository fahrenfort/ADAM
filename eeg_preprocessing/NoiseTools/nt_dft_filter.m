function y=nt_dft_filter(x,transfer,N)
%y=nt_dft_filter(x,transfer,N) - apply filter using DFT
%
%  y: result
%  
%  x: data (time*channels)
%  transfer: transfer function (size=N/2)
%  N: number of samples to use in DFT (default: next power of 2)
%
% Data are zero-padded to N.

if nargin<3||isempty(N); N=2.^nextpow2(size(x,1)); end
if nargin<2; error('!'); end

if mod(N,2)==1; error('N must be multiple of 2'); end % for convenience

% transfer function - if pair of numbers interpret as bandpass
transfer=transfer(:);
if numel(transfer)==2;
    lo=1+round(N/2*transfer(1)); 
    hi=1+round(N/2*transfer(2));
    if min(lo,hi)<1; error('band limits too low'); end
    if max(lo,hi)>N/2; error('band limits too high'); end
    transfer=zeros(N/2,1);
    transfer(lo:hi)=1;
end

[m,n,o]=size(x);
x=reshape(x,m,n*o); % all trials and channels in parallel

% pad with zeros
x=[x;zeros(N-m,n*o)];

% DFT, apply transfer function, inverse DFT
y=fft(x);
transfer=repmat(transfer,1,n*o);
size(transfer);
y=y.*[transfer;flipud(transfer)];
y=real(ifft(y));

% trim, reshape to original geometry
y=y(1:m,:);
y=reshape(y,[m,n,o]);


