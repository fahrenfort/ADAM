% demo code for MCCA

disp('Effect of the eigenspectrum of the noise');

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
Each data matrix is obtained by mixing 1 sinusoidal target source (the same
for all data matrices) and 10 independent white Gaussian noise sources
(different for each data matrix). 

The noise is full-rank so target and noise can only be separated approximately.
The quality of this approximation depends on the eigenspectrum of the noise.  
%}

nsamples=100000;
nchans=10;
nsets=10;

factors=[1 .6 .4 .3]; % determines slope of noise eigenspectrum
SNRs=10.^-(0:10); % SNRs at which to test MCCA

% create in advance the noise and target sources and mixing matrices
target=sin(2*pi*(1:nsamples)/nsamples)'; 
noise_sources={};
noise_mix={};
target_mix={};
for iSet=1:nsets
    nn=randn(nsamples,nchans);
    noise_sources{iSet}=nt_normcol(nt_pca(nn)); % whiten 
    [noise_mix{iSet}, ~]=qr(randn(nchans)); % random rotation to preserve eigenspectrum    
    target_mix{iSet}=nt_normcol(randn(1,nchans));
end

figure(1); clf;


for iFactor=1:numel(factors);
    factor=factors(iFactor);  % determines slope of noise eigenspectrum
    for iSNR=1:numel(SNRs)
        SNR=SNRs(iSNR);
        dataset=zeros(nsamples,nchans,nsets);
        for iSet=1:nsets
            noise=noise_sources{iSet};
            noise=noise*diag(factor.^(0:(nchans-1))); % impose decreasing profile
            nn=noise*noise_mix{iSet};
            nn=nn/sqrt(mean(nn(:).^2)); % set overall power to 1
            tt=target*target_mix{iSet};
            tt=target/sqrt(mean(tt(:).^2)); % set overall power to 1
            dataset(:,:,iSet)=nn+tt*sqrt(SNR);
        end

        % apply MCCA
        x=dataset(:,:); % concatenate channelwise
        C=x'*x;
        [A,score,AA]=nt_mcca(C,nchans);
        z=x*A;

        r(iSNR)=corr(target,z(:,1));
        rr(:,iSNR)=corr(target,z);
    end
    subplot 121; semilogy(mean(noise.^2), '.-'); 
    ylim([10^-10, 2]);
    title('noise eigenspectrum');
    xlabel('PC'); ylabel('variance'); 
    hold on;
    subplot 122; plot(abs(r), '.-');
    ylim([-.1 1.1]);
    ylabel('correlation'); xlabel('target SNR');
    title('correlation of first SC with target');
    set(gca,'xticklabel',SNRs);
    hold on
    drawnow
        
end
% Fig1: MCCA is more successful at extracting weak targets 
% if the eigenspectrum of the noise is decreasing



