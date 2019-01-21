% demo code for MCCA

disp('hybrid between CorCA and MCCA');

%{
EEG data in response to repeated presentation of a 1 kHz tone for 25 subjects.

Data are organized as a 4D matrix with dimensions time, channel, trial, subject.

Data are preprocessed prior to epoching. Bad channels are detected automatically and removed, 
and data are downsampled to 128 Hz and filtered or detrended (several versions).

%}

clear
fname='../../DATA/MCCA_EXAMPLE_DATA/localiser_hpf_1Hz.mat';
if exist(fname, 'file')
    load (fname); % loads x sr preStim postStim 
else
    error('Download example data from http://audition.ens.fr/adc/NoiseTools/DATA/MCCA_EXAMPLE_DATA/');
end

[nsamples,nchans,ntrials,nsets]=size(x);
t=linspace(-preStim/sr, postStim/sr, nsamples);

figure(1); clf
figure(2); clf
figure(3); clf
figure(4); clf

% for iSet=1:nsets;
%     idx=nt_find_outlier_trials(x(:,:,:,iSet),3);
%     x(:,:,setdiff(1:ntrials,idx),iSet)=0;
%     idx=nt_find_outlier_trials(x(:,:,:,iSet),2);
%     x(:,:,setdiff(1:ntrials,idx),iSet)=0;
% end

for iScramble = 0%[0 1]

    if iScramble; nrepeats=10; else nrepeats=1; end
    
    for iRepeat=1:nrepeats
        
        if iScramble
            y=nt_phase_scramble(x(:,:));
            y=reshape(y,[nsamples,nchans,ntrials,nsets]);
        else y=x; end
            
        % remove mean 
        for iSet=1:nsets
            y(:,:,:,iSet)=nt_demean(y(:,:,:,iSet)); % baseline correction works less well
        end

        % initial dimensionality reduction of each subject's data with PCA
        NKEEP=50; 
        y1=zeros(nsamples,NKEEP,ntrials,nsets);
        for iSet=1:nsets
            a=nt_pca(y(:,:,:,iSet)); 
            a=a(:,1:NKEEP,:); % discard PCs beyond NKEEP
            y1(:,:,:,iSet)=a;
        end

        tmp=permute(y1,[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP,nsets]);
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP);
%         figure(1); clf; %hold on;
%         plot(score(:),'.-', 'linewidth', 2-iScramble); 
        
        NKEEP2=200;
        y2=zeros(size(y1));
        for iSet=1:nsets
            b=pinv(AA{iSet});
            y2(:,:,:,iSet)=nt_mmat(y1(:,:,:,iSet),AA{iSet}(:,1:NKEEP2)*b(1:NKEEP2,:));
        end
        
        NKEEP3=45;
        y3=zeros(nsamples,NKEEP,ntrials,nsets);
        for iSet=1:nsets
            C0=nt_cov(y1(:,:,:,iSet));
            C1=nt_cov(y2(:,:,:,iSet));
            [todss,pwr0,pwr1]=nt_dss0(C0,C1);
            fromdss=pinv(todss);
            y3(:,:,:,iSet)=nt_mmat(y2(:,:,:,iSet),todss(:,1:NKEEP3)*fromdss(1:NKEEP3,:));
        end
        
        NKEEP4=40;
        y4=zeros(nsamples,NKEEP4,ntrials,nsets);
        y5=zeros(nsamples,NKEEP4,ntrials,nsets);
        for iSet=1:nsets
            a=nt_pca(y1(:,:,:,iSet));
            y4(:,:,:,iSet)=a(:,1:NKEEP4,:);
            a=nt_pca(y3(:,:,:,iSet));
            y5(:,:,:,iSet)=a(:,1:NKEEP4,:);
        end
            
        
        
        tmp=permute(y4,[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP4,nsets]);
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP4);
        figure(1); hold on;
        plot(score(1:40),'.-', 'linewidth', 2-iScramble); 

        z=tmp*A;
        z=nt_fold(z,nsamples);
        figure(2); clf;
        nt_bsplot(z(:,1,:)); 
        drawnow

        tmp=permute(y5,[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP4,nsets]);
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP4);
        figure(3); hold on;
        plot(score(1:40),'.-', 'linewidth', 2-iScramble); 
        
        z=tmp*A;
        z=nt_fold(z,nsamples);
        figure(4); clf;
        nt_bsplot(z(:,1,:));
        drawnow
        
        
      end
end

