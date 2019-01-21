%{

Extract component based on reproducible induced activity

%}

reset(RandStream.getGlobalStream); % to get reproducible signals

sr=1000;
nsamples=1000;
ntrials=100;
nchans=10;
CF1=4;
CF2=2; 
FSIZE=12;

% sinusoidal pulse 
target=sin(2*pi*0.5*(1:nsamples)'/sr).^2;
target=repmat(target,[1,1,ntrials]) .* sin( 2*pi*(CF1*repmat((1:nsamples)',[1,1,ntrials])/sr + repmat(rand(1,1,ntrials),[nsamples,1,1])));
target=nt_normcol(target);

% noise
NNOISE=9;
noise=randn(nsamples,NNOISE,ntrials);
noise=nt_mmat(noise,randn(NNOISE,nchans));
noise=nt_normcol(noise);

% data
SNR=0.0001;
x=noise+ SNR * nt_mmat(target,randn(1,nchans)) ;

x=nt_demean(x);
x=nt_normcol(x);

[squares,quads,D]=nt_qca(x,[],[],10);


figure(1); clf;
set(gca,'fontsize',12)
set(gcf,'color',[1 1 1])
set(gcf, 'position', [667   368   700   400])

subplot 231
plot(0.9*squeeze(target(:,1,1:5)/max(target(:))), 'k'); set(gca,'ytick',[]); title('target (5 trials)', 'fontsize',14);
subplot 232
plot(squeeze(x(:,1,:))); set(gca,'ytick',[]); title('mixture', 'fontsize',14); xlabel('samples', 'fontsize', 14);
subplot 233
plot(0.9*squeeze(squares(:,1,1:5)/max(squares(:))), 'k');  set(gca,'ytick',[]); title('recovered', 'fontsize',14);

clear x; clear noise


%return
P=mfilename('fullpath');
[PATHSTR,NAME,EXT] = fileparts(P);

% noise is meg data
load([PATHSTR,'/../DATA/meg']);
meg=nt_unfold(meg);
meg=nt_fold(meg(1:76000,:),1000);
meg=nt_demean2(meg);
sr=1000; % Hz

[idx,d]=nt_find_outlier_trials2(meg,1.5);% plot(d);
meg=meg(:,:,idx);
[idx,d]=nt_find_outlier_trials2(meg,1.5); % plot(d);
meg=meg(:,:,idx);
nsamples=size(meg,1); nchans=size(meg,2); ntrials=size(meg,3);


% add 'target'
if 0
    target=sin(2*pi*0.5*(1:nsamples)'/sr).^2;
    target=repmat(target,[1,1,ntrials]).*sin(2*pi*(20*repmat((1:nsamples)',[1,1,ntrials])/sr+repmat(rand(1,1,ntrials),[nsamples,1,1])));
else
    % target is pulses with random latency
    PWIDTH=60; PRANGE=400;
    B=sin(2*pi*(1:PWIDTH)/(2*PWIDTH));
    target=zeros(nsamples,1,ntrials);
    mixmatrix=randn(1,nchans);
    for k=1:ntrials
        target(200+round(rand*PRANGE),:,k)=1;
        target(:,:,k)=filter(B,1,target(:,:,k));
    end
end

SNR=0.001;
mixmatrix=randn(1,nchans);
[dummy,maxidx]=max(abs(mixmatrix));
target=nt_mmat(target,mixmatrix);
meg=nt_normcol(meg)+SNR*target;


y=nt_normcol(nt_pca(nt_normcol(meg)));
NKEEP=28;
yy=zeros(size(meg,1),NKEEP*(NKEEP+1),size(meg,3));
ii=1;
for k=1:NKEEP;
    for j=1:k;
        yy(:,ii,:)=meg(:,k,:).*meg(:,j,:);
        ii=ii+1;
    end
end

yy=nt_demean2(yy);
SMOOTH=100;
yy=filter(ones(1,SMOOTH),1,yy);
yy=yy(SMOOTH:end,:,:);

[todss,pwr0,pwr1]=nt_dss1(yy,[],[],0);
z=nt_mmat(yy,todss);


if 0
    c0=nt_cov(meg(SMOOTH:end,:,:));
    c1=nt_cov(meg(SMOOTH:end,:,:).*repmat(z(:,1,:),[1,nchans,1]));
else
    %c0=nt_cov([y1,y2]);
    c0=nt_cov(meg(SMOOTH:end,:,:).*repmat(max(0,z(:,1,:)),[1,nchans,1]));
    c1=nt_cov(meg(SMOOTH:end,:,:).*repmat(min(0,z(:,1,:)),[1,nchans,1]));
end
   
[todss2,pwr0,pwr1]=nt_dss0(c0,c1);
z2=nt_mmat(meg,todss2);

t=(0:size(target,1)-1)/sr;
subplot 234; 
plot(t, squeeze(target(:,1,1:10))/max(max(target(:,1,:))), 'k'); title('target (10 trials)', 'fontsize', 14); set(gca,'ytick',[]); 
xlim([t(1) t(end)]); ylim([-1.2 1.2]);
subplot 235; 
plot(t, squeeze(meg(:,maxidx,:)));  title('mixture', 'fontsize', 14); xlabel('s', 'fontsize',14); set(gca,'ytick',[]); xlim([t(1) t(end)]);
subplot 236; 
plot(t, -squeeze(z2(:,end,1:10))/max(max(abs(z2(:,end,:)))), 'k');  title('recovered', 'fontsize', 14); set(gca,'ytick',[]); 
xlim([t(1) t(end)]); ylim([-1.2 1.2]);


figure(2); clf
plot(D);
title('quality of fit to a square');
xlabel('component'); ylabel('score');



