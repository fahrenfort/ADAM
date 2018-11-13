% demo code for MCCA

disp('Sensitivity to target mismatch between sets');
rng(1);

clear
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
%}

%{
Each data matrix is obtained by mixing a sinusoidal target 
with 9 independent Gaussian noise sources (different for each data matrix). 

The noise power is either white or in 1/f^2 (obtained by
applying cumsum() to white noise --> random walk).

The target, nominally the same for all data matrices, is mixed with a
controlled amount of Gaussian noise to simulate between-set mismatch.

%}

nsamples=10000;
nchans=20;
Ks=[ 0.01 .1 1]; % amplitude of data matrix-specific mismatch

f=10; 
target=sin(2*pi*f*(1:nsamples)'/nsamples);
target=nt_normcol(target);
target(1:round(nsamples/f*3))=0; target=flipud(target);
target(1:round(nsamples/f*3))=0; % --> sinusoidal pulse

figure(1); clf;
nt_banner('background noise is white');
for iK=1:numel(Ks)
    colored=0; % white background
    [T,recovered]=doit(target,Ks(iK),colored,nchans);
    subplot (3,2,1+2*(iK-1));
    plot(T); set(gca,'ytick',[]); title(['target (K= ',num2str(Ks(iK)), ')']);
    subplot (3,2,2+2*(iK-1)); 
    plot(recovered); set(gca,'ytick',[]); title('SC 1');
end
% Fig1: if the 'noise' is white, the target is recovered (and even enhanced) even if it is not
% identical across sets.

figure(2); clf;
nt_banner('background is 1/f^2');
for iK=1:numel(Ks)
    colored=1; % 1/f^2 background
    [T,recovered]=doit(target,Ks(iK),colored,nchans);
    subplot (3,2,1+2*(iK-1));
    plot(T); set(gca,'ytick',[]); title(['target (K= ',num2str(Ks(iK)), ')']);
    subplot (3,2,2+2*(iK-1)); 
    plot(recovered); set(gca,'ytick',[]); title('SC 1');
end
% Fig2: if the 'noise' is colored, MCCA fails to recover the target

Ks=[.01 .02 .05 .1  .2 .5 1];
Ns=[5 10 20 40 80];
a=[]; aa=[];
for iK=1:numel(Ks);
    for iN=1:numel(Ns);
        [T,recovered]=doit(target,Ks(iK),0,Ns(iN));
        a(iK,iN)=abs(nt_normcol(T)'*nt_normcol(recovered))/nsamples;
        [T,recovered]=doit(target,Ks(iK),1,Ns(iN));
        aa(iK,iN)=abs(nt_normcol(T)'*nt_normcol(recovered))/nsamples;
    end
end
figure(3); clf
subplot 121; nt_imagescc(a); 
c=colorbar('location','northoutside');
set(get(c,'label'),'string', 'correlation with target'); set(c,'limits',[0 1])
set(gca,'xticklabel',Ns,'yticklabel',Ks); title('background white');
xlabel('number of channels'); ylabel('amplitude of data matrix-specific mismatch');
subplot 122; nt_imagescc(aa);
c=colorbar('location','northoutside');
set(get(c,'label'),'string', 'correlation with target'); set(c,'limits',[0 1])
set(gca,'xticklabel',Ns,'yticklabel',Ks); title('background colored');
xlabel('number of channels'); ylabel('amplitude of data matrix-specific mismatch');


% simulate recovering a target
function [T,recovered]=doit(target,K,colored,nchans)
    % target: nominal target pattern, K: amplitude of data matrix-specific
    % mismatch, colored: if non zero background noise is 1/f^2, else white
    nsets=10;
    SNR=10^0;
    nsamples=size(target,1);
    dataset=zeros(nsamples,nchans,nsets);
    for iSet=1:nsets
        if colored
            noise=nt_normcol(cumsum(nt_demean(randn(nsamples,nchans-1))));
        else
            noise=nt_normcol((nt_demean(randn(nsamples,nchans-1))));
        end
        n=noise*randn(nchans-1,nchans);
        T=target + K*randn(size(target)); % set-specific pattern
        tt=T*randn(1,nchans); 
        dataset(:,:,iSet)=n/rms(n(:))+sqrt(SNR)*tt/rms(tt(:));
    end

    % apply MCCA
    x=dataset(:,:);
    x=nt_demean(x);
    C=x'*x;
    [A,score,AA]=nt_mcca(C,nchans);
    z=x*A;
    recovered=z(:,1);
    % T: target used for last data matrix, recovered: SC 1
end
