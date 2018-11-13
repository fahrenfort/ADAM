% demo code for MCCA

disp('set-specific denoising using CCs');

clear
close all
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 9 data matrices, each with 10 channels.
Each data matrix is obtained by mixing 3 different sinusoidal target sources
and 7 independent white Gaussian noise sources. 
One target source is mixed into one set of 3 matrices, the second into a different 
set of 3, and so-on. The noise sources are different for each data matrix.
Targets and noise are mixed via a random mixing matrix, different for each
set.
%}

nsamples=100000;
nchans=10;
nsets=9;
SNR=10^-8; % target SNR in power, same for all data matrices

% which source (number) goes to which channel (rank)
source_map=[1 1 1 2 2 2 3 3 3];

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
freqs=(1:3)';
sources=sin(2*pi*freqs*(1:nsamples)/nsamples)';
sources=nt_normcol(sources);
for iSet=1:nsets
    source=sources(:,source_map(:,iSet));
    dataset(:,:,iSet)= [source*sqrt(SNR), nt_normcol(randn(nsamples,nchans-size(source,2)))] * randn(10);
end

% apply MCCA
x=dataset(:,:); % concatenate channelwise
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);

% SCs
z=x*A; % SCs

figure(1); clf
subplot 331; plot(sources(:,1)); title('target 1'); set(gca,'ytick',[]);
subplot 334; plot(sources(:,2)); title('target 2'); set(gca,'ytick',[]);
subplot 337; plot(sources(:,3)); title('target 3'); set(gca,'ytick',[]);
subplot 232
plot(score, '.-');
title('variance per SC');
ylabel('variance'); xlabel('SC');
subplot 133
plot(z(:,1:3)); title('SCs 1 to 3'); xlabel('sample'); set(gca,'ytick',[]);
subplot 235
nt_imagescc(nt_normcol(sources)'*nt_normcol(z(:,1:20))/nsamples);
xlabel('SC'); ylabel('source');
title('source-SC correlation');
% Fig1: SCs span the subspace containing the sources (no 1-to-1 mapping)


% CCs
for iSet=1:nsets
    zz{iSet}=dataset(:,:,iSet)*AA{iSet}; 
end

figure(2); clf
for iTarget=1:3
    subplot (3,3,3*iTarget-2);
    % correlation of this target with each of the SCs of one of the sets
    % that contains it
    iSet=3*iTarget-1;
    nt_imagescc(nt_normcol(sources)'*nt_normcol(zz{iSet}(:,1:20))/nsamples);
    set(gca,'xtick',[]); ylabel('target'); 
    title(['correlation with CC{',num2str(iSet),'}']); xlabel('CC')
    subplot (3,3,3*iTarget-1);
    % first 3 SCs for this set
    plot(zz{3*iTarget-1}(:,1:3));
    set(gca,'xtick',[],'ytick',[]);
    title(['CC{',num2str(iSet),'}(1:3)']);
    subplot (3,3,3*iTarget);
    % recovered target
    fromccs=pinv(AA{iSet});
    denoise_matrix=AA{iSet}(:,1:3)*fromccs(1:3,:);
    denoised=dataset(:,:,iSet)*denoise_matrix;
    topcs=nt_pca0(denoised);
    recovered=denoised*topcs(:,1);
    plot(recovered);  
    set(gca,'xtick',[],'ytick',[]);
    title(['recovered target ', num2str(iTarget)]);
end
%Fig2: individual sources can be recovered from the CCs (denoising)

