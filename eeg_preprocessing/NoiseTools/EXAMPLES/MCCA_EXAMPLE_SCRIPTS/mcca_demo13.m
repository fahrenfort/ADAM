% demo code for MCCA

disp('Phase shifted surrogate data, real target');

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
SNR=10^-1;

f=10; 
target=sin(2*pi*f*(1:nsamples)'/nsamples);
target=nt_normcol(target);
target(1:round(nsamples/f*3))=0; target=flipud(target);
target(1:round(nsamples/f*3))=0; 

% original data
figure(1); clf
nt_banner('original data');
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
subplot 121;
plot(score, '.-'); title('1/f^2 noise, SNR=10^-2');
ylabel('variance'); xlabel('SC'); title('variance profile')
subplot 122; 
plot(z(:,1)); title('SC 1')


% phase-shifted surrogate data
figure(2); clf
nt_banner('surrogate, random phase-shift');
for iRepeat=1:10
    for iSet=1:nsets
        dataset(:,:,iSet)=nt_phase_scramble(dataset(:,:,iSet));
    end
    % apply MCCA
    x=dataset(:,:); % concatenate channelwise
    x=nt_demean(x);
    C=x'*x;
    [A,score,AA]=nt_mcca(C,nchans);
    z=x*A;
    subplot 121;
    plot(score, '.-'); title('1/f^2 noise, SNR=10^-2');
    ylabel('variance'); xlabel('SC'); title('variance profile')
    hold on
    subplot 122; 
    plot(z(:,1)); title('SC 1')
    hold on
end
%Fig2: phase randomization destroys the temporal pattern found by MCCA (SC1) in the
%original data, suggesting that that pattern was real.


