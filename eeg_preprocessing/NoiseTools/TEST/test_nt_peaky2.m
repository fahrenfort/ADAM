% Test of nt_peaky_DA on the EEG data

%clear
%close all

%dname='C:\Users\Trevor\Documents\Data\';
dname='/data/meg/arzounian/DA_131218/';

% fnames={
%     'DA_131218_p11\DA_131218_p11_001_spontan_open'
%     'DA_131218_p11\DA_131218_p11_002_spontan_closed'
%     'DA_131218_p11\DA_131218_p11_003_spontan_closed_off'
%     };
fnames={
    'DA_131218_p14_001_spontan_open'
    'DA_131218_002_spontan_closed'
    'DA_131218_003_spontan_closed_off'
    };

start=[65 40 120];
stop=[925 760 1900];

iFile=3;
fname=fnames{iFile};

disp('read data...');
h=sopen([dname,fname]);
sr=h.SampleRate;
x=sread(h);
sclose(h);
disp('done')

% exclude last channels (unused), first and last chunks (artifacted)
x=x(start(iFile)*sr:stop(iFile)*sr,1:37);    

% high-pass
HPF=15;
x=nt_demean(x,1:round(sr/HPF));
[B,A]=butter(2,HPF/(sr/2),'high');
x=filter(B,A,x);
x=x(round(sr/HPF)+1:end,:);


% down-sample
DSR=sr/256;
x=nt_dsample(x,DSR);
sr=sr/DSR;


window=30*sr;   % effective integration window size(samples)
T=10*sr; % Time resolution (samples)
nSmooth=floor(window/T);

tic
x0=x;
x=nt_pca(x);
[tocomps,ii]=nt_peaky([],x,T,nSmooth);
[~,idx]=sort(ii);
%tocomps=tocomps(:,idx);

zz=nt_mmat(x,tocomps);
zz=nt_normcol(zz);

figure(1);
imagescc((filter(ones(round(10*sr),1),1,abs(zz)).^0.5)'); title('all')
title('Reconstruction')
toc

figure(2); clf
subplot 121; plot(kurtosis(zz), '.-'); title('kurtosis'); xlabel('component')
subplot 122; imagescc(nt_cov(zz)); title('covariance matrix'); xlabel('component'); ylabel('component');

figure(3); clf
imagescc(nt_xcov(nt_normcol(x0),nt_normcol(zz))); 
xlabel('component'); ylabel('electrode'); title('electrode/component covariance');

figure(4); clf;
nt_spect_plot2(nt_dsample(zz,2),1024,[],[],sr/2)
title('spectra of components');


