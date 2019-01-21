% demo code for MCCA

disp('MCCA for targets shared 2 by 2');

clear
close all
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
Each data matrix is obtained by mixing 20 different sinusoidal target sources
and 6 independent white Gaussian noise sources. Each target source is mixed
into 2 channels.  The noise sources are different for each data matrix.
Targets and noise are mixed via a random mixing matrix, different for each
set.
%}

nsamples=100000;
nchans=10;
nsets=10;
SNR=10^-4; % target SNR in power, same for all data matrices

% which source goes to which channel (easier to hardwire than calculate)
source_map=[...
    1 1 2 3 4;
    2 5 5 6 7;
    3 6 8 8 9;
    4 7 9 10 10];
source_map =[source_map, source_map+10]; 

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
freqs=(1:20)';
sources=sin(2*pi*freqs*(1:nsamples)/nsamples)';
sources=nt_normcol(sources);
for iSet=1:nsets
    source=sources(:,source_map(:,iSet));
    dataset(:,:,iSet)= [source*sqrt(SNR), nt_normcol(randn(nsamples,nsets-size(source,2)))] * randn(10);
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
% Fig1: there are more shared components than channels in the data

figure(2); clf;
subplot 211;
plot(z(:,1:20)); title('SCs 1 to 20'); xlabel('sample');
subplot 212;
plot(z(:,21)); title('SC 21'); xlabel('sample');
% Fig2: MCCA correctly isolates a subspace spanning the shared sources

figure(3); clf;
nt_imagescc(nt_normcol(sources)'*nt_normcol(z(:,1:20))/nsamples);
xlabel('SC'); ylabel('source');
title('correlation between sources and SCs');
h=colorbar('limits',[-1 1]); set(get(h,'label'),'string','correlation')
% Fig3: SCs span the subspace containing the sources (no 1-to-1 mapping)


