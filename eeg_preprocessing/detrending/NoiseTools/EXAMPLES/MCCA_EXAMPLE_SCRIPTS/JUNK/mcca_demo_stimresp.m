%rng(1)
clc; clear; close all
add_uncorrelated_noise = 0.000; % <------ tweak this
nsamples=10000;
nchans=5;
nsets=5;
SNR=10^-10; % target SNR in power, same for all data matrices


%{
The synthetic dataset consists of 5 data matrices, each with 10 channels.
Each data matrix is obtained by mixing 1 Gaussian target source (the same
for all data matrices) and 9 independent white Gaussian noise sources
(different for each data matrix).
Mixing matrices are random and different for each data matrix.

A stimulus-response analysis is used to quantify the relation between target
and synthetic data. For illustration purposes, the data is here split into a
training set (5000 samples) a validation set (2500 samples) and a test set
(2500 samples). It is assumed that there is an instantaneous coupling
between target and synthetic responses, so no lags are considered here.
Ridge regression is use to train filters that can reconstruct the target
based on the training- and validation set or alternatively predict the data 
based on the target input. In all cases, we evaluate how MCCA denoising
affects the stimulus-response mappings
%}

% dataset is structured as a 3D matrix
dataset=zeros(nsamples,nchans,nsets);
target=randn(nsamples,1);
target=nt_normcol(target); % normalize for convenience
for iSet=1:nsets
    noise=randn(nsamples,nchans-1); % noise sources
    noise=nt_normcol(noise); % normalize for convenience
    mix_noise=nt_normcol(randn(nchans-1,nchans));
    mix_target=nt_normcol(randn(1,nchans))*sqrt(SNR);
    dataset(:,:,iSet)=target*mix_target+noise*mix_noise;
    
    if add_uncorrelated_noise
        dataset(:,:,iSet) = dataset(:,:,iSet)+ nt_normcol(randn(nsamples,nchans))*sqrt(10^(-10));
    end
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
[A,score,AA] = nt_mcca(C,nchans);

% for later use, we will also train decoders on the mcca rotated data
z_dataset = zeros(nsamples,nchans*nsets,nsets);
for iSet=1:nsets
    z_dataset(:,:,iSet) = dataset(:,:,iSet)*AA{iSet};
end



% we use ridge regression to estimate the stimulus-response relation.
nlambda = 50;
lambda = logspace(-13,13,nlambda);
ridgeestimator = @(xx,yy,ll) (xx'*xx + ll*eye(size(xx,2)))\(xx'*yy);
% initialize correlation matrices
cvalts = zeros(nsets,2);  cvaltr = zeros(nlambda,nsets); cvaltr_mcca = zeros(nlambda,nchans*nsets,nsets);
fcvalts = zeros(nsets,2); fcvaltr = zeros(nlambda,nchans,nsets); fcvaltr_mcca = zeros(nlambda,nchans*nsets,nsets);
maxfun = @(x)(max(x(:))); 


% *************************************************************************
% *************************************************************************

% 1) train decoders that can extract the target based on raw per-subject data  
for iSet=1:nsets
    
    % use the training set and validation set to estimate an optimal ridge
    % parameter
    for l = 1 : nlambda
        W = ridgeestimator(dataset(ixtr,:,iSet),target(ixtr),lambda(l));
        cvaltr(l,iSet) =  corr(dataset(ixvl,:,iSet)*W,target(ixvl),'Type','Pearson');
    end
    
    % the ridge parameter that gives best validation-set score is selected
    [~,bestlambda] = max(cvaltr(:,iSet));
    W = ridgeestimator(dataset([ixtr ixvl],:,iSet),target([ixtr ixvl]),lambda(bestlambda));
    
    % evaluate how well the "optimal" ridge filter "reconstructs" the target
    cvalts(iSet,1) = corr(target(ixts),dataset(ixts,:,iSet)*W,'Type','Pearson');
end


% 2) train decoders that can extract the target based on "mcca denoised" per-subject data  
for iSet=1:nsets
    
    % although not strictly necessary, we do a grid search over all
    % lambda values and over all subspace values (i.e. number of mcca
    % components to retain). we seek the most parsimonious model that
    % gives the highest correlation coefficients over the validation set
    for n = 1 : nchans*nsets
        for l = 1 : nlambda
            W = ridgeestimator(z_dataset(ixtr,1:n,iSet),target(ixtr),lambda(l));
            cvaltr_mcca(l,n,iSet) =  corr(z_dataset(ixvl,1:n,iSet)*W,target(ixvl),'Type','Pearson');
        end
    end
    [~,idmax] = maxfun(cvaltr_mcca(:,:,iSet));
    [lm,nm]=ind2sub([nlambda,nchans*nsets],idmax);
    
    
    % train the decoder
    W = ridgeestimator(z_dataset([ixtr ixvl],1:nm,iSet),target([ixtr ixvl]),lambda(lm));
    
    
    % evaluate performance
    cvalts(iSet,2) = corr(z_dataset(ixts,1:nm,iSet)*W,target(ixts),'Type','Pearson');
end



% 3) train encoders that can predict the raw per-subject data from target
for iSet=1:nsets

    % loop over all channels
    for ch = 1 : nchans
        for l = 1 : nlambda
            Wf = ridgeestimator(target(ixtr),dataset(ixtr,ch,iSet),lambda(l));
            cvaltr(l,ch,iSet) =  corr(dataset(ixvl,ch,iSet),target(ixvl)*Wf,'Type','Pearson');
        end
        
        % the ridge parameter that gives best validation-set score is selected
        [~,bestlambda] = max(cvaltr(:,ch,iSet));
        Wf = ridgeestimator(target([ixtr ixvl]),dataset([ixtr ixvl],ch,iSet),lambda(bestlambda));
        
        % evaluate how well the "optimal" ridge filter predicts raw data
        fcvalts(iSet,ch,1) = corr(target(ixts)*Wf,dataset(ixts,ch,iSet),'Type','Pearson');
    end
end


% 4) train encoders that can predict the "denoised" per-subject data from target
for iSet=1:nsets
    
     for ch = 1 : nchans*nsets
        for l = 1 : nlambda
            Wf = ridgeestimator(target(ixtr),z_dataset(ixtr,ch,iSet),lambda(l));
            fcvaltr_mcca(l,ch,iSet) =  corr(z_dataset(ixvl,ch,iSet),target(ixvl)*Wf,'Type','Pearson');
        end
        [~,bestlambda] = max(fcvaltr_mcca(:,ch,iSet));
        Wf = ridgeestimator(target([ixtr ixvl]),z_dataset([ixtr ixvl],ch,iSet),lambda(bestlambda));
        
        % evaluate how well the "optimal" ridge filter predicts "denoised" data
       fcvalts(iSet,ch,2) = corr(target(ixts)*Wf,z_dataset(ixts,ch,iSet),'Type','Pearson');
    end
end



figure('Position',[100 100 500 800]); clf
bar(cvalts)
ylabel('Correlation between reconstruction and target')
xlabel('Subject')
hleg = legend('Without MCCA','With MCCA')
hleg.Box = 'off';
hleg.Location = 'southoutside';
title({'Reconstructing the target based','on per-subject stimulus-response models'},'FontWeight','Normal')



figure('Position',[800 100 400 200]); clf
subplot 121
imagesc(fcvalts(:,1:nchans,1),[-1 1])
xlabel('Channels'); ylabel('Subjects')
title({'Predicting','channel responses'},'FontWeight','normal')
subplot 122
imagesc(fcvalts(:,1:nchans,2),[-1 1])
xlabel('MCCA components'); ylabel('Subjects')
title({'Predicting','component responses'},'FontWeight','normal')