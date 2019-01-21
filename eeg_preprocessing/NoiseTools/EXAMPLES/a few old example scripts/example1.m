% Find the linear combinations of multichannel data that
% maximize repeatability over trials.  Data are time*channel*trials.
%
% Uses nt_dss0().

clear;
disp(mfilename);
help(mfilename)

% create synthetic data
nsamples=100*3;
nchans=30;
ntrials=100;
noise_dim=20; % dimensionality of noise
source=[zeros(nsamples/3,1);sin(2*pi*(1:nsamples/3)/(nsamples/3))';zeros(nsamples/3,1)]; 
s=source*randn(1,nchans);
s=repmat(s,[1,1,100]); % evoked
SNR=0.1;
noise=nt_mmat(randn(nsamples,noise_dim,ntrials), randn(noise_dim,nchans));
data=noise/rms(noise(:))+SNR*s/nt_rms(s(:));

% apply DSS to clean them
c0=nt_cov(data);
c1=nt_cov(mean(data,3));
[todss,pwr0,pwr1]=nt_dss0(c0,c1);
z=nt_mmat(data,todss);


% plot results
figure(1); clf
subplot 131; 
plot(source); title('source'); 
subplot 132;
plot(mean(data,3)); title('data');
subplot 133;
nt_bsplot(z(:,1,:)); title('recovered'); 

