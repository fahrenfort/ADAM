% demo code for MCCA

disp('Competition between real and artifactual components');

clear
% check for NoiseTools
if 2 ~= exist('nt_mcca')
    error('Download NoiseTools from http://audition.ens.fr/adc/NoiseTools/ and put on path');
end

%{
The synthetic dataset consists of 10 data matrices, each with 10 channels.
%}

%{
Each data matrix is obtained by mixing a sinusoidal target (same for each data matrix) 
with 9 independent Gaussian noise sources
(different for each data matrix). The noise power is in 1/f^2 (obtained by
applying cumsum() to white noise --> random walk).
%}

nsamples=100000;
nchans=10;
nsets=10;
SNR=10^-4;

f=10; 
target=sin(2*pi*f*(1:nsamples)'/nsamples);
target=nt_normcol(target);
target(1:round(nsamples/f*3))=0; target=flipud(target);
target(1:round(nsamples/f*3))=0; 

figure(1); clf

% white noise
dataset=zeros(nsamples,nchans,nsets);
for iSet=1:nsets
    noise=(nt_demean(randn(nsamples,nchans)));
    n=noise*randn(nchans,nchans);
    t=target*randn(1,nchans);
    dataset(:,:,iSet)=n/rms(n(:))+sqrt(SNR)*t/rms(t(:));
end
% apply MCCA
x=dataset(:,:); % concatenate channelwise
x=nt_demean(x);
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;
subplot 331;
plot(score, '.-');  title('white noise')
ylabel('variance'); xlabel('SC');
subplot 332; 
plot(z(:,1)); title('SC 1')
subplot 333;
plot(abs(corr(target,z)), '.-'); ylim([-.05 1.05])
xlabel('SC'); ylabel('abs(r)'); title('correlation with target');

% 1/f^2 noise, SNR=10^-4
dataset=zeros(nsamples,nchans,nsets);
for iSet=1:nsets
    noise=cumsum(nt_demean(randn(nsamples,nchans)));
    n=noise*randn(nchans,nchans);
    t=target*randn(1,nchans);
    dataset(:,:,iSet)=n/rms(n(:))+sqrt(SNR)*t/rms(t(:));
end
% apply MCCA
x=dataset(:,:); % concatenate channelwise
x=nt_demean(x);
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;
subplot 334;
plot(score, '.-'); title('1/f^2 noise, SNR=10^-4');
ylabel('variance'); xlabel('SC');
subplot 335; 
plot(z(:,1)); title('SC 1')
subplot 336;
plot(abs(corr(target,z)), '.-'); ylim([-.05 1.05])
xlabel('SC'); ylabel('abs(r)'); title('correlation with target');

% 1/f^2 noise, SNR=10^-2
SNR=10^-2;
dataset=zeros(nsamples,nchans,nsets);
for iSet=1:nsets
    noise=cumsum(nt_demean(randn(nsamples,nchans)));
    n=noise*randn(nchans,nchans);
    t=target*randn(1,nchans);
    dataset(:,:,iSet)=n/rms(n(:))+sqrt(SNR)*t/rms(t(:));
end
% apply MCCA
x=dataset(:,:); % concatenate channelwise
x=nt_demean(x);
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;
subplot 337;
plot(score, '.-'); title('1/f^2 noise, SNR=10^-2');
ylabel('variance'); xlabel('SC');
subplot 338; 
plot(z(:,1)); title('SC 1')
subplot 339;
plot(abs(corr(target,z)), '.-'); ylim([-.05 1.05])
xlabel('SC'); ylabel('abs(r)'); title('correlation with target');



