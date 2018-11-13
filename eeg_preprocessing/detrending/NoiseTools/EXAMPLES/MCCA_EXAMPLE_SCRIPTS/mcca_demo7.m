% demo code for MCCA

disp('MCCA is biased towards low-frequency components');

clear
close all
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
%}

%{
Each data matrix is obtained by mixing 10 independent Gaussian noise sources
(different for each data matrix). The noise power is in 1/f^2 (obtained by
applying cumsum() to white noise --> random walk).

No target.
%}

nsamples=100000;
nchans=10;
nsets=10;

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
for iSet=1:nsets
    noise=randn(nsamples,nchans); % 
    noise=cumsum(nt_demean(noise));
    dataset(:,:,iSet)=noise;
end

figure(1); clf

% apply MCCA
x=dataset(:,:); % concatenate channelwise
x=nt_demean(x);
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;
plot(score, '.-');
ylabel('variance'); xlabel('SC');
% Fig1: The first few SCs have high scores, despite the lack of any shared
% target. This is due to spurious correlations between low-pass noise
% components.

figure(2); clf
subplot 221
sr=1000; % Hz - arbitrary sampling rate
dsr=1; nt_spect_plot(nt_resample(x,1,dsr),256,[],[],sr/dsr);
title('power spectrum of raw data');
subplot 223
dsr=64; nt_spect_plot2(nt_resample(z,1,dsr),256,[],[],sr/dsr)
ylabel('SC');
title('power spectrum for each SC');
scs=[1 2 10 100];
for iSC=1:4;
    subplot(4,2,iSC*2);
    plot(z(:,scs(iSC))); 
    set(gca,'ytick',[],'xtick',[]);
    title(['SC ', num2str(scs(iSC))]);
end
set(gca,'xtick',[0,nsamples]); xlabel('samples')
% Fig2: SCs are biased to low frequencies.  The first SC has the lowest
% frequency content, followed by the second, etc.  Each SC has a different
% spectrum, despite the fact that all data channels had the same (low pass)
% spectrum.

