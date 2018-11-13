function [c0,c1]=nt_bias_fft(x,freq,nfft)
%[c0,c1]=nt_bias_fft(x,freq,nfft) - covariance with and w/o filter bias
%
% x: data 
% freq: row vector of normalized frequencies to keep (wrt sr)
% nfft: fft size
%
% The filter has zeros at all frequencies except those immediately inferior
% or superior to values in vector freq.
% 
% If freq has two rows, keep frequencies between corresponding values on
% first and second row.
%
% NoiseTools

if max(freq(:))>0.5; error('frequencies should be <= 0.5'); end
if nfft>size(x,1); error('nfft too large'); end

filt=zeros(nfft/2+1,1);

if size(freq,1)==1
    for k=1:size(freq,2)
        idx=round(freq(1,k)*nfft+0.5);
        filt(idx)=1;
    end
elseif size(freq,1)==2
    for k=1:size(freq,2)
        idx=round(freq(1,k)*nfft+0.5) : round(freq(2,k)*nfft+0.5);
        filt(idx)=1;
    end
else
    error('freq should have one or two rows');
end

filt=[filt;flipud(filt(2:end-1))];

%plot(filt); pause

w=hanning(nfft);

%plot(filt); return

[m,n,o]=size(x);
c0=nt_cov(x);
c1=zeros(size(c0));
for j=1:o
    nframes=ceil((m-nfft/2)/(nfft/2));
    for k=1:nframes
        idx=(k-1)*nfft/2;
        idx=min(idx,m-nfft);
        z=x(idx+1:idx+nfft,:,j);
        Z=fft(nt_vecmult(z,w));
        Z=nt_vecmult(Z,filt);
        c1=c1+real(Z'*Z);
    end
end

%[todss,fromdss,ratio,pwr]=dss0(c0,c1);     
