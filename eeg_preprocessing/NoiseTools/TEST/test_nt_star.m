clear

%{
Test 1:
N channels,  1 sinusoidal target, N-3 noise sources, temporally local
artifacts on each channel.
%}

nsamples=1000;
f=2;
target=sin((1:nsamples)'/nsamples*2*pi*f);
nchans=10;
noise=randn(nsamples,nchans-3);

SNR=sqrt(1);
x0=nt_normcol(noise*randn(size(noise,2),nchans))+ SNR * target*randn(1,nchans);
x0=nt_demean(x0);
artifact=zeros(size(x0));
for k=1:nchans
    artifact( (k-1)*100 + (1:20) , k)=1;
end
x=x0+20*artifact;

[y,w]=nt_star(x,2); 
figure(1); clf; subplot 311; plot(x); subplot 312;  plot(y); title(['SNR= ' ,num2str(SNR)]);
subplot 313; plot(nt_demean(y-x0));

SNR=sqrt(10^-7);
x0=nt_normcol(noise*randn(size(noise,2),nchans))+ SNR * target*randn(1,nchans);
x0=nt_demean(x0);
artifact=zeros(size(x0));
for k=1:nchans
    artifact( (k-1)*100 + (1:20) , k)=1;
end
x=x0+20*artifact;

%x=[x,x(:,1)]; % make data rank deficient

[y,w]=nt_star(x,2); 

figure(2); clf; subplot 311; plot(x);title(['SNR= ' ,num2str(SNR)]);
subplot 312;  plot(y); title('artifact removed');
subplot 313; plot(nt_demean(y-x0)); title('error after artifact removal')


x=nt_demean(x);
y=nt_demean(y);

c0=nt_cov(x);
c1=nt_cov(nt_detrend(x,2,[],'sinusoids'));
todss=nt_dss0(c0,c1);
z1=nt_normcol(nt_mmat(x,todss));
c0=nt_cov(y);
c1=nt_cov(nt_detrend(y,2,[],'sinusoids'));
todss=nt_dss0(c0,c1);
z2=nt_normcol(nt_mmat(y,todss));

figure(3); clf; set(gcf, 'color',[1 1 1], 'position', [300   500   520 600])

subplot 511
plot(target)
set(gca,'fontsize', 14, 'xtick',[], 'ytick', []);
ylim([-1.1 1.1]); title('target')

subplot 512
plot(x)
set(gca,'fontsize', 14, 'xtick',[], 'ytick', []);
mn=min(x(:)); mx=max(x(:));
ylim([mn-(mx-mn)*0.1, mx+(mx-mn)*0.1])
title('mixed with noise and artifacts')

subplot 513
plot(y)
set(gca,'fontsize', 14, 'xtick',[], 'ytick', []);
title('STAR')

subplot 515
plot(z2(:,end));
mn=min(z2(:,end)); mx=max(z2(:,end));
set(gca,'fontsize', 14, 'ytick', []); xlabel('samples');
ylim([mn-(mx-mn)*0.1, mx+(mx-mn)*0.1])
%text(800,mx-(mx-mn)*0.1, 'SNR=10^-8','fontsize',14)
title('STAR + JD')

subplot 514
plot(z1(:,end));
set(gca,'fontsize', 14, 'xtick',[], 'ytick', []);
mn=min(z1(:,end)); mx=max(z1(:,end));
ylim([mn-(mx-mn)*0.1, mx+(mx-mn)*0.1])
title('JD')




