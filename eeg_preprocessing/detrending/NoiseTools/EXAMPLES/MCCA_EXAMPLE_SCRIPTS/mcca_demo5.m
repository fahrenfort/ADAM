% mcca_demo5 - demo code for MCCA

disp('MCCA with rank deficient data.')
disp('Applying dimensionality reduction before MCCA avoids numerical issues.');

clear
close all
%rng('default'); % reset to same state to get deterministic result
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
%}

%{
Each data matrix is obtained by mixing 1 sinusoidal target source (the same
for all data matrices) and 8 independent white Gaussian noise sources
(different for each data matrix). 
Each data matrix is thus rank-deficient.
%}

nsamples=100000;
nchans=5;
nsets=10;
SNR=10^-20; % target SNR in power, same for all data matrices

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
target=sin(2*pi*(1:nsamples)/nsamples)'; % 1 target source
target=nt_normcol(target); % normalize for convenience
for iSet=1:nsets
    NNOISE=nchans-2;
    noise=cumsum(randn(nsamples,NNOISE)); % 8 noise sources
    noise=nt_normcol(noise); % normalize for convenience
    forward_noise=nt_normcol(randn(NNOISE,nchans));
    forward_target=nt_normcol(randn(1,nchans))*sqrt(SNR);
    dataset(:,:,iSet)=noise*forward_noise + target*forward_target;
end

figure(1); clf

% apply MCCA
x=dataset(:,:); % concatenate channelwise
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z1=x*A;
subplot 121
plot(score, '.-');
title('Without dimensionalitry reduction');
ylabel('variance'); xlabel('SC');

% reduce dimensionality of each data matrix, then apply MCCA
NDISCARD=1;
dataset_reduced=zeros(nsamples,nchans-NDISCARD,nsets);
for iSet=1:nsets
    a=nt_pca(dataset(:,:,iSet));
    dataset_reduced(:,:,iSet)=a(:,1:nchans-NDISCARD);
end
x=dataset_reduced(:,:); % concatenate channelwise
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans-NDISCARD);
z2=x*A;
subplot 122
plot(score, '.-');
title('With dimensionality reduction');
ylabel('variance'); xlabel('SC');
% Fig1: rank deficient data causes numerical problems (left), reflected by
% spurious variance scores. Reducing dimensionality before MCCA cures them
% (right).

figure(2); clf
subplot 121; plot(z1(:,1)); title('SC1, without dim. reduction');
subplot 122; plot(z2(:,1)); title('SC1, with dim. reduction');

