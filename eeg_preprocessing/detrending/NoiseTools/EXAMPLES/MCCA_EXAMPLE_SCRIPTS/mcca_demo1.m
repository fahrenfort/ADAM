% demo code for MCCA

disp('Compare publicly available implementations and simple Matlab code');

clear
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end
% check for Lucas Parra's mcca
if 2 ~= exist('corrca')
    error('Download mcca from https://www.parralab.org/corrca/ and put in path');
else
    [FILEPATH,NAME,EXT] = fileparts(which('corrca'));
    if 2~=exist([FILEPATH,filesep,'mcca'])
        error('Download mcca from https://www.parralab.org/corrca/ and put in path');
    end
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
Each data matrix is obtained by mixing 1 sinusoidal target source (the same
for all data matrices) and 9 independent white Gaussian noise sources
(different for each data matrix).
Mixing matrices are random and different for each data matrix. 
%}

nsamples=100000;
nchans=10;
nsets=10;
% target SNR in power, same for all data matrices:
%SNR=10^-4; % usually works for all
SNR=10^-8; % can fail for mcca.m
%SNR=10^-14; % usually fails for mcca.m

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

if 0
    % plot correlation of each channel with target to show that the target can't be 
    % seen in the raw data (correlation values are consistent with spurious correlations)
    figure(2); clf;
    subplot 121;
    target=nt_normcol(target); % unit norm
    for iSet=1:nsets
        x=dataset(:,:,iSet); 
        x=nt_normcol(x); % unit norm
        r=x'*target/nsamples; % correlation with target
        plot(r, '.-'); 
        hold on;
    end
    title('correlation raw data with target'); xlabel('channel'); ylabel('correlation');
    legend('dataset 1', 'dataset 2', 'etc.'); legend boxoff

    % plot best-correlated channel across all datasets
    x=dataset(:,:); % concatenate channelwise
    r=x'*target;
    [~,idx]=max(abs(r));
    subplot 222
    plot(target); 
    title('target'); xlabel('samples');
    subplot 224
    plot(x(:,idx));
    title('best-correlated channel/dataset'); xlabel('samples');
    % Fig1: target is not evident in raw data (poor SNR)
end

%{
Apply MCCA using several implementations.
%}

figure(1); clf

% NoiseTools
x=dataset(:,:); % concatenate channelwise
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;
subplot 131
plot(z(:,1));
title('nt_mcca.m','interpreter','none'); xlabel('samples');
legend('first component'); % first SC (summary component)
legend boxoff

% Lucas Parra's mcca (from corrca code)
x=dataset(:,:); % concatenate channelwise
d=repmat(nchans,1,nsets);
%[V,rho,A,rhotest,Amean]=mcca(x,d); !!! fails
[V,rho]=mcca(x,d);
z=x*V;
subplot 132
plot(z(:,1));
title('mcca.m'); xlabel('samples');
legend('first component');
legend boxoff


% simple matlab code
y=zeros(size(dataset));
for iSet=1:nsets
    y(:,:,iSet)=nt_normcol(nt_pca(dataset(:,:,iSet))); % whiten 
end
y=y(:,:); % concatenate over channels
z=nt_pca(y);
subplot 133
plot(z(:,1));
title('matlab'); xlabel('samples');
legend('first component');
legend boxoff


%{
Summary: MCCA can be implemented in several ways.
Implementations may differ in robustness, and the results may have different scaling. 
See documentation of each implementation for further details and features.
%}



    







    
