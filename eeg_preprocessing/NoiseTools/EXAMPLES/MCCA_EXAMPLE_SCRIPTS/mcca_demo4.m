% demo code for MCCA

disp('MCCA for multiple target sources');

clear
close all
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
Each data matrix is obtained by mixing 2 sinusoidal target sources (the same
for all data matrices) and 8 independent white Gaussian noise sources
(different for each data matrix).
Mixing matrices are random and different for each data matrix. 
%}

nsamples=100000;
nchans=10;
nsets=10;
SNR=10^-2; % target SNR in power, same for all data matrices

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
freqs=[1,2]';
target=sin(2*pi*freqs*(1:nsamples)/nsamples)'; % target sources
target=nt_normcol(target); % normalize for convenience
for iSet=1:nsets
    noise=randn(nsamples,nchans-2); % noise sources
    noise=nt_normcol(noise); % normalize for convenience
    forward_noise=nt_normcol(randn(nchans-2,nchans));
    forward_target=nt_normcol(randn(2,nchans))*sqrt(SNR);
    dataset(:,:,iSet)=noise*forward_noise + target*forward_target;
end

% apply MCCA
x=dataset(:,:); % concatenate channelwise
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;

figure(1); clf
plot(score, '.-');
title('variance per SC');
ylabel('variance'); xlabel('SC');
% Fig1: variance profile

figure(2); clf;
subplot 121
plot(target); 
title('targets'); xlabel('samples');
subplot 122
plot(z(:,1:2));
title('first two SCs'); xlabel('samples');
% Fig2: each SC is a linear combination of targets


figure(3); clf;
plot(nt_normcol(target)'*nt_normcol(z(:,1:2))/nsamples, '.-');
xlabel('target source'); ylabel('correlation');
legend('SC1', 'SC2'); legend boxoff
title('correlation between sources and SCs');
xlim([0 3]); ylim([-1 1]);
set(gca,'xtick', [1 2])
% Fig3: SCs span the subspace containing targets (no 1-to-1 mapping)

