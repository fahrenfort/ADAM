% demo code for MCCA

disp('Random phase-shifted surrogate data');

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

% original data
figure(1); clf
nt_banner('original data')
subplot 121
x=dataset(:,:); % concatenate channelwise
x=nt_demean(x);
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;
plot(score, '.-');
ylabel('variance'); xlabel('SC');
title('variance profile');
subplot 122
plot(z(:,1)); title('SC1');

% surrogate phase-shifted data
figure(2); 
nt_banner('surrogate, random phase-shift');
for iRepeat=1:10;
    for iSet=1:nsets
        dataset(:,:,iSet)=nt_phase_scramble(dataset(:,:,iSet));
    end
    x=dataset(:,:); % concatenate channelwise
    x=nt_demean(x);
    C=x'*x;
    [A,score,AA]=nt_mcca(C,nchans);
    z=x*A;
    subplot 121
    plot(score, '.-');
    ylabel('variance'); xlabel('SC'); 
    title('variance profile');
    hold on
    subplot 122
    plot(z(:,1));
    title('SC 1');
    hold on
end
% Fig2: random phase shifting of each dataset should destroy any trace of a
% real shared component.  Nonetheless, each analysis shows much the same
% variance profile, and a clearly-shaped first SC component.  

