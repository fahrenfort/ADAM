% demo code for MCCA

disp('Show MCCA transform matrices'); 

clear
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end
rng(1);

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
Each data matrix is obtained by mixing 1 sinusoidal target source (the same
for all data matrices) and 9 independent white Gaussian noise sources
(different for each data matrix).
Mixing matrices are random and different for each data matrix. 
%}

nsamples=100000;
nchans=10;
nsets=5;
SNR=10^-12; % target SNR in power, same for all data matrices

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
target=sin(2*pi*(1:nsamples)/nsamples)'; % target source
target=nt_normcol(target); % normalize for convenience
for iSet=1:nsets
    noise=randn(nsamples,nchans-1); % noise sources
    noise=nt_normcol(noise); % normalize for convenience
    mix_noise=nt_normcol(randn(nchans-1,nchans));
    mix_target=nt_normcol(randn(1,nchans))*sqrt(SNR);
    dataset(:,:,iSet)=target*mix_target+noise*mix_noise;
end

% apply MCCA
x=dataset(:,:); % concatenate channelwise
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;

%{
output arguments are:
A: transform matrix from concatenated data to summary components (SC)
score: variance of each SC
AA: array of transform matrices from each data matrix to canonical components (CC)
%}

figure(1); clf;
subplot 121
imagesc(A);
set(gca,'clim', [-max(abs(A(:))),max(abs(A(:)))]);
title('SC transform matrix A');
xlabel('SC'); ylabel('channel (concatenated data)');
for iSet=1:nsets
    subplot(nsets,2,iSet*2);
    imagesc(AA{iSet});
    set(gca,'clim', [-max(abs(A(:))),max(abs(A(:)))]);
    set(gca,'xtick',[], 'ytick', [1,nchans]);
    title(['CC transform matrix AA{',num2str(iSet),'}'], 'interpreter', 'none');
    ylabel('channel')
end
set(gca,'xtick',0:10:100);
xlabel('CC');
% Fig2: CC transform matrices AA{} are slices of the SC transform matrix A 


% plot SCs and CCs
zz={};
for iSet=1:nsets
    zz{iSet}=dataset(:,:,iSet)*AA{iSet};
end

figure(2); clf
plot (z(:,1), 'r','linewidth', 2); hold on
for iSet=1:nsets
    plot(zz{iSet}(:,1));
end
title('each SC is the sum of CCs'); 
xlabel('sample');
legend('SC (:,1)', 'CC\{1\} (:,1)', 'CC\{2\} (:,1)', 'etc.'); legend boxoff
% Fig3: each SC is the sum of CCs across data matrices

