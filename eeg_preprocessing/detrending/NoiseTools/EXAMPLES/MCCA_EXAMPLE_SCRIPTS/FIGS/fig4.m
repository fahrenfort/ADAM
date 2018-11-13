% figures for mcca paper

nchans=10;
nnoise=9;
nsamples=10000;
nsets=10;
SNR=10^-10; % signal to noise ratio (in power)

nn=randn(nsamples,nnoise*nsets);
nn=reshape(nn,[nsamples,nnoise,nsets]);

% sinusoidal target
target=sin(2*pi*[0:nsamples-1]/nsamples)';
%target=randn(size(target));

% synthesize all signals
xx=zeros(nsamples,nchans,nsets);
for k=1:nsets
    xx(:,:,k) = sqrt(SNR) * target*randn(1,nchans)+ nn(:,:,k) * randn(nnoise,nchans);
end
xx0=xx;
xx0=reshape(xx0,[nsamples, nchans*nsets]);

if 0
    xx=reshape(xx,[nsamples, nchans*nsets]);
    c=nt_cov(xx)/size(xx,1);
    A=nt_mcca(c,nchans);
    y=xx*A;
else
    % whiten spatially each dataset
    for k=1:nsets
        xx(:,:,k)=nt_normcol(nt_pca(xx(:,:,k)));
    end
    xx=reshape(xx,[nsamples, nchans*nsets]);
    y=nt_pca(xx);
end
    
    
% xx=reshape(xx,[nsamples, nchans*nsets]);
% nt_imagescc(xx); 
% return

    


figure(1); clf; set(gcf,'color', [1 1 1], 'position',[440   585   600   200])
axes('position', [.1 .2 .16 .70]);
dstime=1:30:size(target,1);
plot(dstime,target(dstime,:));
set(gca,'fontsize',14, 'ytick',[]); title('target')
ylim(1.5*max(target)*[-1 1])
xlabel('sample')

axes('position', [.32 .2 .16 .70])
plot(xx0(:,1:10));
set(gca,'fontsize',14, 'ytick',[]); title('target + noise')
ylim(max(xx0(:))*[-1 1])
xlabel('sample')

axes('position', [.56 .2 .16 .70])
p=mean(y.^2);
plot(mean(y.^2), '.-k','markersize',22);
set(gca,'fontsize',14, 'ytick',[0 5 10]); xlabel('SC'); ylabel ('variance');
ylim([0 12])


axes('position', [.78 .2 .16 .70])
plot(y(:,1));
set(gca,'fontsize',14, 'ytick',[]); title('recovered')
ylim(1.5*max(y(:,1))*[-1 1])
xlabel('sample')

save Fig_A_data p



return


figure(2); clf; nt_imagescc(xx)
figure(3); clf; 
[~,idx]=max(abs(nt_xcov(target,xx))); 
plot(nt_normcol([xx(:,idx),y(:,1)])); legend('best pc','mcca')