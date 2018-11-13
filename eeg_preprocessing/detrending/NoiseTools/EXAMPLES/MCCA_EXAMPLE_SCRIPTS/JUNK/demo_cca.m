
%rng(1)
clc; clear; close all
nsamples=1000;
nchans=5;
nsets=10;
SNR=10^-3; % target SNR in power, same for all data matrices
lagoi = 0:32;
dataset=zeros(nsamples,nchans,nsets);


target=randn(nsamples,1);

% predefine an encoding filter
w_encoder = nt_normcol([normpdf(0:32,8,2)*0.2-normpdf(0:32,15,2) + 0.3*normpdf(0:32,22,2)]');

% shift target and normalize
target_shifted = [nt_multishift(target,lagoi); zeros(lagoi(end),length(lagoi))];
target_source = nt_normcol(target_shifted*w_encoder);



for iSet=1:nsets
    noise=randn(nsamples,nchans-1); % noise sources
    noise=nt_normcol(noise); % normalize for convenience
    mix_noise=nt_normcol(randn(nchans-1,nchans));
    mix_target=nt_normcol(randn(1,nchans))*sqrt(SNR);
    
    

    % this is quick and very dirty code! A BETTER EXAMPLE SHOULD BE
    % DESIGNED!
    % here, we assume that the target source spuriously mix out to the
    % channel responses with random mixing matrices in 200 sample long
    % blocks. these mixing matrices are random within-subject
    % and across subjects. 
    wdur = 200;
    win = 1 : wdur; mix_distractor =[];
    while win(end)<=nsamples
        distractor = target_shifted(win,:)*w_encoder;
        if randi(2)==1
            mix_distractor = [mix_distractor; nt_normcol(distractor)*nt_normcol(randn(1,nchans))]; 
        else
            mix_distractor = [mix_distractor; zeros(wdur,nchans)];
        end
        win = win + wdur;
    end
    
    dataset(:,:,iSet)=target_source*mix_target + nt_normcol(mix_distractor+ noise*mix_noise);
    
    
    
end

% we split the dataset into a training set, a validation set and a test set
ixtr = 1 : nsamples/2;
ixvl = ixtr(end)+1 : (ixtr(end)+1 + nsamples/4);
ixts = ixvl(end)+1 : nsamples;

% standardize the dataset based on the mean and standard deviation
% estimated from the training set
for iSet=1:nsets
    dataset(:,:,iSet) = bsxfun(@rdivide,bsxfun(@minus,dataset(:,:,iSet),mean(dataset(ixtr,:,iSet))),...
        std(dataset(ixtr,:,iSet)));
end

% train the mcca matrices over the training set
x = dataset(ixtr,:,:);  x = x(:,:);% concatenate channelwise
C = x'*x;
[A,score,AA] = nt_mcca(C,nchans); % <- currently we DON'T optimize prewhitening K 
z_dataset = zeros(nsamples,nchans*nsets,nsets);
for iSet=1:nsets
    z_dataset(:,:,iSet) = dataset(:,:,iSet)*AA{iSet};
end



% when should the "last" mcca components be removed????? sometimes they
% contain relevant information
% z_dataset = z_dataset(:,1:size(dataset,2),:); <-------------------



for with_without_mcca = 1 : 2
    if with_without_mcca == 1
        fprintf('\n fitting models to the untransformed data')
    else
        fprintf('\n fitting models to the MCCA denoised data')
    end
    for iSet=1:nsets
        
        if with_without_mcca ==1
            yy = dataset(:,:,iSet);
        else
            yy = z_dataset(:,:,iSet);
        end
        
        fprintf('.')
        
        
        xtr = target_shifted(ixtr,:);
        ytr = yy(ixtr,:);
        xvl = target_shifted(ixvl,:);
        yvl = yy(ixvl,:);
        xts = target_shifted(ixts,:);
        yts = yy(ixts,:);
        
        clear xtopcs ytopcs cc_first_component w_encoder_recovered w_decoder_recovered
        [xtopcs]=nt_pca0(xtr);
        [ytopcs]=nt_pca0(ytr);
        for nx = 1 : size(xtopcs)
            for ny = 1 : size(ytopcs)
                [A,B]=nt_cca(xtr*xtopcs(:,1:nx),ytr*ytopcs(:,1:ny));
                
                w_encoder_recovered = xtopcs(:,1:nx)*A;
                w_decoder_recovered = ytopcs(:,1:ny)*B;
                cc_first_component(nx,ny)=corr2(xvl*w_encoder_recovered(:,1),yvl*w_decoder_recovered(:,1)).^2;
            end
        end
        [~,best_nx_ny] = max(cc_first_component(:));
        [nxm,nym]=ind2sub([size(cc_first_component,1),size(cc_first_component,2)],best_nx_ny);
        
        
        % now, based on the optimal set of Kx and Ky values, fit the CCA
        % model again
        [A,B]=nt_cca([xtr; xvl]*xtopcs(:,1:nxm),[ytr; yvl]*ytopcs(:,1:nym));

        % the below filters should now be optimal
        w_encoder_recovered = xtopcs(:,1:nxm)*A;
        w_decoder_recovered = ytopcs(:,1:nym)*B;
        
        
        w_encoder_recovered_all(:,iSet,with_without_mcca) = -sign(w_encoder_recovered(15,1))*w_encoder_recovered(:,1);
        cc_test_set(iSet,with_without_mcca)=corr2(xts*w_encoder_recovered(:,1),yts*w_decoder_recovered(:,1)).^2;
    end
end


disp(cc_test_set)

figure('Position',[100 100 600 700]); clf
subplot(2,3,[1 2 4 5 ]);
bar(mean(cc_test_set))
hold on
plot(cc_test_set','Color',[0.9 0.9 0.9])
set(gca,'xtick',[1 2]); set(gca,'xticklabel',{'Without MCCA','With MCCA'})
ylabel('Squared correlation coefficient between first CC pair')


subplot(2,3,[3]);
plot(zscore(w_encoder),'.-','Color',[0.6 0.6 0.6])
hold on
plot(mean(zscore(w_encoder_recovered_all(:,:,1)),2),'-k')
title({'Recovered filter','(without MCCA)'})
set(gca,'xtick',[]); set(gca,'ytick',[])
hleg = legend('Target','1st component');
hleg.Box = 'off';
hleg.Box = 'off';
hleg.Location = 'southoutside';

subplot(2,3,[6]);
plot(zscore(w_encoder),'.-','Color',[0.6 0.6 0.6])
hold on
plot(mean(zscore(w_encoder_recovered_all(:,:,2)),2),'-k')
title({'Recovered filter','(with MCCA)'})
hleg = legend('Target','1st component');
hleg.Box = 'off';
hleg.Location = 'southoutside';
set(gca,'xtick',[]); set(gca,'ytick',[])








