% demo code for MCCA

disp('Benefit of dimensionality reduction (synthetic data)');
rng(1);

clear
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

nsamples=10000;
nchans=40;

f=10; 
target=sin(2*pi*f*(1:nsamples)'/nsamples);
target=nt_normcol(target);
target(1:round(nsamples/f*3))=0; target=flipud(target);
target(1:round(nsamples/f*3))=0; % --> sinusoidal pulse


%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
%}

%{
Each data matrix is obtained by mixing a sinusoidal target 
with 9 independent Gaussian noise sources (different for each data matrix). 

The noise power is in 1/f^2 (obtained by applying cumsum() to white noise --> random walk).
The variance profile of the  noise is shaped (high-order components are
weaker).

The target, nominally the same for all data matrices, is mixed with
Gaussian noise to simulate between-set differences. The amplitude of this
data matrix-specific mismatch is varied as a parameter.

Dimensionality reduction is applied to each data matrix (PCA then select).
Number of PCs retained is varied as a parameter.

%}


Ks=[.05 .07 .1 .14 .2 .28];
nPCs=[6 8 10 12 14 16 18 20 22 24 26 30];
a=[]; aa=[];


SNR=1; factor=0.7;
for iK=1:numel(Ks);
    for iNkeep=1:numel(nPCs);
        [T,recovered]=doit(target,Ks(iK),nchans,nPCs(iNkeep),SNR,factor);
        aa(iK,iNkeep)=abs(nt_normcol(T)'*nt_normcol(recovered))/nsamples;
    end
end
figure(1); clf
nt_imagescc(aa);
c=colorbar('location','eastoutside');
set(get(c,'label'),'string', 'correlation with target'); set(c,'limits',[0 1])
set(gca,'xtick',1:numel(nPCs),'xticklabel',nPCs,'ytick',1:numel(Ks),'yticklabel',Ks); title('SNR=1; factor=0.7');
xlabel('number of PCs retained'); ylabel('amplitude of data matrix-specific mismatch');

SNR=10^-2; factor=0.8;
for iK=1:numel(Ks);
    for iNkeep=1:numel(nPCs);
        [T,recovered]=doit(target,Ks(iK),nchans,nPCs(iNkeep),SNR,factor);
        aa(iK,iNkeep)=abs(nt_normcol(T)'*nt_normcol(recovered))/nsamples;
    end
end
figure(2); clf
nt_imagescc(aa);
c=colorbar('location','eastoutside');
set(get(c,'label'),'string', 'correlation with target'); set(c,'limits',[0 1])
set(gca,'xtick',1:numel(nPCs),'xticklabel',nPCs,'ytick',1:numel(Ks),'yticklabel',Ks); title('SNR=10^{-2}; factor=0.8');
xlabel('number of PCs retained'); ylabel('amplitude of data matrix-specific mismatch');

SNR=10^-2; factor=0.9;
for iK=1:numel(Ks);
    for iNkeep=1:numel(nPCs);
        [T,recovered]=doit(target,Ks(iK),nchans,nPCs(iNkeep),SNR,factor);
        aa(iK,iNkeep)=abs(nt_normcol(T)'*nt_normcol(recovered))/nsamples;
    end
end
figure(3); clf
nt_imagescc(aa);
c=colorbar('location','eastoutside');
set(get(c,'label'),'string', 'correlation with target'); set(c,'limits',[0 1])
set(gca,'xtick',1:numel(nPCs),'xticklabel',nPCs,'ytick',1:numel(Ks),'yticklabel',Ks); title('SNR=10^{-2}; factor=0.9');
xlabel('number of PCs retained'); ylabel('amplitude of data matrix-specific mismatch');


% simulate recovering a target
function [T,recovered]=doit(target,K,nchans,nPCs,SNR,factor)
    % target: nominal target pattern, K: amplitude of data matrix-specific
    % mismatch, nchans: number of channels, nPCs: number of PCs to retain
    % factor: determines steepness of noise variance profile
    rng(1);
    nsets=10;
    nsamples=size(target,1);
    dataset=zeros(nsamples,nchans,nsets);
    for iSet=1:nsets
        noise=nt_normcol(cumsum(nt_demean(randn(nsamples,nchans-1))));
        for k=1:size(noise,2);
            noise(:,k)=factor.^k * noise(:,k);  
        end
        n=noise*randn(nchans-1,nchans);
        T=target + K*randn(size(target)); % set-specific pattern
        tt=T*randn(1,nchans); 
        dataset(:,:,iSet)=n/rms(n(:))+sqrt(SNR)*tt/rms(tt(:));
    end
    
    dataset2=zeros(nsamples, nPCs, nsets);
    for iSet=1:nsets;
        a=nt_pca(dataset(:,:,iSet));
        dataset2(:,:,iSet)=a(:,1:nPCs);
    end
    dataset=dataset2;

    % apply MCCA
    x=dataset(:,:);
    x=nt_demean(x);
    C=x'*x;
    [A,score,AA]=nt_mcca(C,nPCs);
    z=x*A;
    recovered=z(:,1);
    % T: target used for last data matrix, recovered: SC 1
end
