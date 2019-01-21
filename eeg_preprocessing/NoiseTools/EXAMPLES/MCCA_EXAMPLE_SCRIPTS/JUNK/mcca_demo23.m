% demo code for MCCA

disp('Real EEG data: MCCA and DSS in tandem');

%{
EEG data in response to repeated presentation of a 1 kHz tone for 25 subjects.

Data are organized as a 4D matrix with dimensions time, channel, trial, subject.

Data are preprocessed prior to epoching. Bad channels are detected automatically and removed, 
and data are downsampled to 128 Hz and filtered or detrended (several versions).

Analysis alternates between MCCA and DSS.  Each is used to denoise and
reduce dimensionality to make the other more effective.

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
figure(5); clf

for iSet=1:nsets;
    idx=nt_find_outlier_trials(x(:,:,:,iSet),3);
    x(:,:,setdiff(1:ntrials,idx),iSet)=0;
    idx=nt_find_outlier_trials(x(:,:,:,iSet),2);
    x(:,:,setdiff(1:ntrials,idx),iSet)=0;
end

for iScramble =[0 1]

    if iScramble; nrepeats=10; else nrepeats=1; end
    
    for iRepeat=1:nrepeats
                
        if iScramble
            y=nt_phase_scramble(x(:,:));
            x=reshape(y,[nsamples,nchans,ntrials,nsets]);
        end

        % remove mean over pre-stimulus segment ('baseline correction')
        for iSet=1:nsets
            x(:,:,:,iSet)=nt_demean(x(:,:,:,iSet),find(t<0));
        end

        % initial dimensionality reduction of each subject's data with PCA
        NKEEP=40; 
        y1=zeros(nsamples,NKEEP,ntrials,nsets);
        for iSet=1:nsets
            a=nt_pca(x(:,:,:,iSet)); 
            a=a(:,1:NKEEP,:);
            y1(:,:,:,iSet)=a;
        end

        y2=y1;
        
        % first MCCA of all data, use it to reduce dimensionality of each subject's data
        tmp=permute(y2,[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP,nsets]);
        disp(size(tmp))
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP);
        
        figure(2); hold on; %clf;
        plot(score(1:20),'.-', 'linewidth', 2-iScramble); 
        NKEEP3=100;
        NKEEP4=30;
        y3=zeros(nsamples,NKEEP4,ntrials,nsets);
        for iSet=1:nsets;
            back=pinv(AA{iSet});
            denoise_matrix=AA{iSet}(:,1:NKEEP3)*back(1:NKEEP3,:);
            a=nt_mmat(y2(:,:,:,iSet), denoise_matrix);
            a=nt_pca(a);
            y3(:,:,:,iSet)=a(:,1:NKEEP4,:);
        end

        % second DSS
        figure(3); hold on; %clf;
        NKEEP5=12; 
        y4=zeros(nsamples,NKEEP5,ntrials,nsets);
        for iSet=1:nsets;
            [todss,pwr0,pwr1]=nt_dss1(y3(:,:,:,iSet));
            y4(:,:,:,iSet)=nt_mmat(y3(:,:,:,iSet),todss(:,1:NKEEP5));
            subplot(5,5,iSet);
            plot(pwr1(1:10)./pwr0(1:10),'.-'); hold on; drawnow
        end

        % second MCCA
        tmp=permute(y4,[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP5,nsets]);
        tmp=mean(y4,3);
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP5);
        
        figure(4); hold on; %clf;
        plot(score(1:20),'.-', 'linewidth', 2-iScramble); 
        
        NKEEP6=10;
        y5=zeros(nsamples, NKEEP6, ntrials, nsets);
        for iSet=1:nsets
            denoise_matrix=AA{iSet}(:,1:NKEEP6);
            y5(:,:,:,iSet)=nt_mmat(y4(:,:,:,iSet),denoise_matrix);
        end
        
        % final DSS
        figure(5); 
        for iSet=1:nsets
            [todss,pwr0,pwr1]=nt_dss1(y5(:,:,:,iSet));
            %y2(:,:,:,iSet)=nt_mmat(y1(:,:,:,iSet),todss(:,1:NKEEP2));
            subplot(5,5,iSet);
            plot(pwr1./pwr0,'.-'); hold on;
            drawnow
        end
    end
end

