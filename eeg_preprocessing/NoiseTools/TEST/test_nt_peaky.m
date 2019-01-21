% simulation with multichannel data

%{
Data are uncorrelated noise sources mixed via a random matrix.  One or
more of these sources are slowly modulated.  The aim is to discover this
structure and/or extract the source.
%}
clear
close all

nsamples=10000;
nchannels=10;

if 1
    x=randn(nsamples,nchannels);
    x=nt_normcol(nt_pca(x)); % decorrelate noise sources
else
    x=zeros(nsamples,nchannels);
    freqs=100+30*(1:nchannels);
    for k=1:nchannels
        x(:,k)=sin(2*pi*freqs(k)*(1:nsamples)'/nsamples);
    end
end

% modulate
SNRs=ones(nchannels,1);
depths=ones(nchannels,1);
for k=1:nchannels 
    x(:,k)=SNRs(k)*x(:,k).*(1+depths(k)*sin(2*pi*(rand+k*((1:nsamples)'/nsamples)))); 
end    
x0=x; % save

figure(1); clf
subplot 311
imagescc(filter(ones(10,1),1,abs(x0))'); title('Sources')

% mix
x=x*randn(nchannels); 
subplot 312
imagescc(filter(ones(10,1),1,abs(x))'); title('Mixture')

figure(2); clf
subplot 211; 
plot(x); title('Mixture')
subplot 212;
p=mean(nt_pca(x).^2);
plot(p/max(p), '.-');set(gca,'yscale','log');
title('PCA spectrum');
xlabel('PC'); ylabel('power')


%%
T=100;  % Time resolution (samples)
W=200;  % Effective integration time
nSmooth=floor(W/T);

[tocomps,score]=nt_peaky([],x,T,nSmooth);

zz=nt_mmat(x,tocomps);
zz=nt_normcol(zz);
figure(1);
subplot 313
imagescc(filter(ones(100,1),1,abs(zz))'); title('all')
title('Reconstruction')
