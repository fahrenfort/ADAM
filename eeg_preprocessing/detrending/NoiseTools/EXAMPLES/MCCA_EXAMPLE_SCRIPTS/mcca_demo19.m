% demo code for MCCA

disp('Sensitivity to imperfect similarity between sets');
rng(1);

clear
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
%}

%{
Each data matrix is obtained by mixing a sinusoidal target 
with 9 independent Gaussian noise sources (different for each data matrix). 

The noise power is either white or in 1/f^2 (obtained by
applying cumsum() to white noise --> random walk).

The target, nominally the same for all data matrices, is mixed with a
controlled amount of Gaussian noise to simulate between-set differences.

%}

nsamples=10000;
nchans=10;
nsets=10;
SNR=10^-4;

f=10; 
target=sin(2*pi*f*(1:nsamples)'/nsamples);
target=nt_normcol(target);
target(1:round(nsamples/f*3))=0; target=flipud(target);
target(1:round(nsamples/f*3))=0; % --> sinusoidal pulse

figure(1); clf;
specificities=[0 0.1 1 3];
for iSpec=1:numel(specificities)
    dataset=zeros(nsamples,nchans,nsets);
    for iSet=1:nsets
        noise=nt_normcol(cumsum(nt_demean(randn(nsamples,nchans-1))));
        n=noise*randn(nchans-1,nchans);
        t=target + specificities(iSpec)*randn(size(target)); % set-specific pattern
        tt=t*randn(1,nchans); 
        dataset(:,:,iSet)=n/rms(n(:))+sqrt(SNR)*tt/rms(tt(:));
    end

    % apply MCCA
    x=dataset(:,:);
    x=nt_demean(x);
    C=x'*x;
    [A,score,AA]=nt_mcca(C,nchans);
    z=x*A;

    subplot (4,3,1+3*(iSpec-1));
    plot(t); set(gca,'ytick',[]); title(['target (spec= ',num2str(specificities(iSpec)), ')']);
    subplot (4,3,2+3*(iSpec-1));
    plot(score, '.-'); title('1/f^2 noise, SNR=10^-2');
    ylabel('variance'); xlabel('SC'); title('variance profile');
    subplot (4,3,3+3*(iSpec-1)); 
    plot(z(:,1)); set(gca,'ytick',[]); title('SC 1');
end

figure(1); clf;
specificities=[0 0.1 1 3];
for iSpec=1:numel(specificities)
    dataset=zeros(nsamples,nchans,nsets);
    for iSet=1:nsets
        noise=nt_normcol(cumsum(nt_demean(randn(nsamples,nchans-1))));
        n=noise*randn(nchans-1,nchans);
        t=target + specificities(iSpec)*randn(size(target)); % set-specific pattern
        tt=t*randn(1,nchans); 
        dataset(:,:,iSet)=n/rms(n(:))+sqrt(SNR)*tt/rms(tt(:));
    end

    % apply MCCA
    x=dataset(:,:);
    x=nt_demean(x);
    C=x'*x;
    [A,score,AA]=nt_mcca(C,nchans);
    z=x*A;

    subplot (4,3,1+3*(iSpec-1));
    plot(t); set(gca,'ytick',[]); title(['target (spec= ',num2str(specificities(iSpec)), ')']);
    subplot (4,3,2+3*(iSpec-1));
    plot(score, '.-'); title('1/f^2 noise, SNR=10^-2');
    ylabel('variance'); xlabel('SC'); title('variance profile');
    subplot (4,3,3+3*(iSpec-1)); 
    plot(z(:,1)); set(gca,'ytick',[]); title('SC 1');
end



