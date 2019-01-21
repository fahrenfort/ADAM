% demo code for MCCA

disp('Estimate mixing matrix from MCCA solution');


clear
close all
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 5 data matrices, each with 10 channels.
Each data matrix is obtained by mixing 1 sinusoidal target source (the same
for all data matrices) and 9 independent white Gaussian noise sources
(different for each data matrix).
Mixing matrices are random and different for each data matrix. 
%}

nsamples=100000;
nchans=10;
nsets=5;
SNR=10^-2; % target SNR in power, same for all data matrices

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
target=sin(2*pi*(1:nsamples)/nsamples)'; % target source
target=nt_normcol(target); % normalize for convenience
for iSet=1:nsets
    noise=randn(nsamples,nchans-1); % noise sources
    noise=nt_normcol(noise); % normalize for convenience
    forward_noise=nt_normcol(randn(nchans-1,nchans));
    forward_target{iSet}=nt_normcol(randn(1,nchans))*sqrt(SNR);
    dataset(:,:,iSet)=target*forward_target{iSet} + noise*forward_noise;
end

% apply MCCA
x=dataset(:,:); % concatenate channelwise
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);

%{
The forward matrix weights for the target in each data matrix can be
estimated as the cross-correlation between the corresponding CC and the raw
data for that dataset.  The target is assumed to be of norm 1.  Sign is
arbitrary.
%}

figure(1); clf;
iSet=1; 
z=dataset(:,:,iSet)*AA{iSet};
a=dataset(:,:,iSet)'*nt_normcol(z(:,1))/nsamples;
a=a.*sign(a'*forward_target{iSet}'); % align signs
plot(forward_target{iSet}, '.-');
hold on;
plot (a, '.-');
legend('mixing weights (data matrix 1)', 'correlation between CC 1 and raw data'); legend boxoff
xlabel('channel');
% Fig3: correlation coefficients between the first CC and the raw data are
% our best estimate of the (hidden) mixing weights for the target source




