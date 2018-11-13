function [s,f,t]=nt_sgram(x,window,noverlap,nfft,sr,flags)
%[s,f,t]=nt_sgram(x,window,noverlap,nfft,sr,flags) - spectrogram
%
% Without output arguments: plot cubic root of power spectrogram.
% With output arguments: return power spectrogram.
%  s: spectrogram matrix (frequency X time)
%  f: array of frequencies
%  t: array of times
%
%  x: data (if multidimensional, power is averaged over higher dimensions)
%  window: size of hamming window in samples, if vector, shape of window
%  noverlap: frame period in samples
%  nfft: number of points in FFT (default 2^nextpow2(window size))
%  sr: sampling rate.
%  flags: 'nomean': remove mean of each window
%
%  NoiseTools
nt_greetings;

if nargin<1; error('!'); end
if nargin<2; window=[]; end
if nargin<3; noverlap=[]; end
if nargin<4; nfft=[]; end
if nargin<5; sr=[]; end
if nargin<6; flags=[]; end

if isempty(window)
    window=2^nextpow2(round(size(x,1)/8));
end
window=window(:);
if numel(window)==1;
    window=hanning(window);
end
if isempty(noverlap)
    noverlap=floor(numel(window)/2);
end
if isempty(nfft)
    nfft=2^nextpow2(numel(window));
end
if isempty(sr)
    sr=1;
end


x=reshape(x,[size(x,1), prod(size(x))/size(x,1)]); % fold all higher dimensions

nframes=floor( (size(x,1)-numel(window))/noverlap )  +  1;
if nframes<1; error('data too short for window size'); end

s=zeros(nfft/2+1,nframes);

for k=1:nframes
    start=(k-1)*noverlap;
    xx=nt_vecmult(x(start+1:start+size(window),:),window); 
    if ~isempty(flags) & strcmp(flags,'demean');
        xx=nt_demean(xx);
    end
    XX=abs(fft(xx,nfft)).^2;
    s(:,k)=sum(XX(1:(nfft/2+1),:),2);
end

f=(0:nfft/2)*sr/nfft;
t= (numel(window)/2 + (0:nframes-1)*noverlap) / sr;



if nargout==0;
    ss=s.^(1/5);
    %ss=ss-repmat(mean(ss,2),1,size(ss,2));

    nt_imagescc(ss);
    ytick=niceticks(sr/2);
    set(gca,'ytick',1+ytick*nfft/sr, 'yticklabel',num2str(ytick'));
    ylabel('Hz');
    xtick=niceticks(size(x,1)/sr);
    set(gca,'xtick',1+xtick/noverlap*sr - nfft/2/noverlap, 'xticklabel',num2str(xtick'));
    xlabel('s');
    set(gca,'tickdir','out')
    s=[];t=[];f=[];
    axis xy
end


