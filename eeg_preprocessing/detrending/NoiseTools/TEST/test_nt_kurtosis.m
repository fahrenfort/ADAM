
% Test nt_bias_kurtosis on synthetic data
clf;
rng(0); % reset random number generator 
nsamples=5000;
nchans=10;

%{
We create a peaky target and mix it with multiple sources of noise via a
random mixing matrix.
%}

target=zeros(nsamples,1);
target(900:1100,:)=1; %target(1001:1100,:)=-1;
noise=randn(nsamples,nchans);
x=noise(:,1:nchans-1)*randn(nchans-1,nchans);
SNR=0.001;
x=x + SNR*target*randn(1,nchans);
x=nt_demean(x);

figure(1); clf
subplot 221; plot (target); title('target'); ylim([0 1.1]);
subplot 222; plot(x); title('mixture');

%{
We use nt_bias_kurtosis to find a transform that reveals components with
high kurtosis.
%}
nIterations=5;
exponent=2;
w=[];
smooth=1;
[todss]=nt_kurtosis(x,nIterations,exponent,w,smooth);

z=nt_mmat(x,todss);
subplot 224
plot(z(:,1)); title('extracted');

%{
Same, with two targets.
%}

target1=zeros(nsamples,1);
target1(900:1000,:)=1; %target(1001:1100,:)=-1;
target2=zeros(nsamples,1);
target2(2001:2200,:)=1; %target(1001:1100,:)=-1;
x=noise(:,1:nchans-2)*randn(nchans-2,nchans);
SNR=0.01;
x=x + SNR*target1*randn(1,nchans) + SNR*target2*randn(1,nchans);
x=nt_demean(x);

figure(2); clf
subplot 221; plot ([target1,target2]); title('targets'); ylim([0 1.1]);
subplot 222; plot(x); title('mixture');

nIterations=5;
exponent=2;
w=[];
smooth=1;
[todss]=nt_kurtosis(x,nIterations,exponent,w,smooth);
z=nt_mmat(x,todss);

subplot 426
plot(z(:,1)); title('first peaky component');

%{
We use the peakiest component to define a weight that deemphasizes the
interval for which that component has high amplitude.  
We then call nt_bias_kurtosis to find remaining peaky components.
%}

w=abs(z(:,1))<0.5*max(abs(z(:,1)));
nIterations=5;
exponent=2;
smooth=1;
[todss]=nt_kurtosis(x,nIterations,exponent,w,smooth);
z=nt_mmat(x,todss);

subplot 428;
plot(z(:,1)); title('second peaky component');

