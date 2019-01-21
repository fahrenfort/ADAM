% demo code for MCCA

disp('MCCA and DSS in tandem');

%{
EEG data in response to repeated presentation of a 1 kHz tone for 25 subjects.

Data are organized as a 4D matrix with dimensions time, channel, trial, subject.

Data are preprocessed prior to epoching. Bad channels are detected automatically and removed, 
and data are downsampled to 128 Hz and filtered or detrended (several versions).

Analysis alternates between MCCA and DSS.  Each is used to denoise and
reduce dimensionality to make the other more effective.

This script should be seen as 'proof of concept'.  No attempt has been made
to optimize the many parameters.

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

for iScramble =[0 1]

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

        % initial MCCA of all data (to plot initial baseline variance profile)
        NKEEP2=40; % same as second to ease comparison
        tmp=permute(y1(:,1:NKEEP2,:,:),[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP2,nsets]);
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP2);
        figure(1); subplot 131
        plot(score(1:20),'.-', 'linewidth', 2-iScramble);
        title('initial MCCA'); xlabel('SC'); ylabel('variance'); ylim([1 4]);  hold on;

        % first DSS of each subject's data
        figure(2); 
        nt_banner('first DSS');
        y2=zeros(nsamples,NKEEP,ntrials,nsets);
        for iSet=1:nsets
            [todss,pwr0,pwr1]=nt_dss1(y1(:,:,:,iSet));
            y2(:,:,:,iSet)=nt_mmat(y1(:,:,:,iSet),todss);
            subplot(5,5,iSet); hold on
            plot(pwr1(1:10)./pwr0(1:10),'.-');
            if iSet==21; xlabel('component'); ylabel('score'); else; set(gca,'xtick',[]); end
            title(iSet);
            drawnow
        end

        % second MCCA of all data, use it to reduce dimensionality of each subject's data
        y2=y2(:,1:NKEEP2,:,:); % discard DSS components beyond NKEEP2
        tmp=permute(y2,[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP2,nsets]);
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP2);        
        figure(1); subplot 132
        plot(score(1:20),'.-', 'linewidth', 2-iScramble); ylim([1 4]); hold on
        title('after first DSS');
        
        % based on second MCCA, denoise each subject
        NKEEP3=100;
        NKEEP4=35;
        y3=zeros(nsamples,NKEEP4,ntrials,nsets);
        for iSet=1:nsets;
            back=pinv(AA{iSet}); % AA{iSet} projects to CCs, back projects back
            D=AA{iSet}(:,1:NKEEP3)*back(1:NKEEP3,:); % discard CCs beyond NKEEP3
            a=nt_mmat(y2(:,:,:,iSet), D);
            a=nt_pca(a); % PCA denoise data to 'concentrate' shared dimensions
            y3(:,:,:,iSet)=a(:,1:NKEEP4,:); % discard PCs behond NKEEP4
        end

        % second DSS
        figure(3);
        nt_banner('second DSS');
        NKEEP5=12; 
        y4=zeros(nsamples,NKEEP5,ntrials,nsets);
        for iSet=1:nsets;
            [todss,pwr0,pwr1]=nt_dss1(y3(:,:,:,iSet));
             subplot(5,5,iSet); hold on;
            plot(pwr1(1:10)./pwr0(1:10),'.-'); 
            if iSet==21; xlabel('component'); ylabel('score'); else; set(gca,'xtick',[]); end
            drawnow;
            title(iSet)
            y4(:,:,:,iSet)=nt_mmat(y3(:,:,:,iSet),todss(:,1:NKEEP5)); % discard DSS components beyond NKEEP5
       end

        % final MCCA
        tmp=permute(y4,[1 3 2 4]);
        tmp=reshape(tmp,[nsamples*ntrials,NKEEP5,nsets]);
        tmp=tmp(:,:);
        C=tmp'*tmp;
        [A,score,AA]=nt_mcca(C,NKEEP5);
        figure(1); subplot 133;
        plot(score(1:20),'.-', 'linewidth', 2-iScramble); ylim([1 4]); hold on;
        title('after second DSS'); legend('intact','phase scrambled'); legend boxoff
        
        % based on final MCCA, denoise each subject
        NKEEP6=10;
        y5=zeros(nsamples, NKEEP6, ntrials, nsets);
        for iSet=1:nsets
            D=AA{iSet}(:,1:NKEEP6); % discard CCs beyond KNEEP6
            y5(:,:,:,iSet)=nt_mmat(y4(:,:,:,iSet),D);
        end
        
        % final DSS
        figure(4); 
        nt_banner('final DSS');
        for iSet=1:nsets
            [todss,pwr0,pwr1]=nt_dss1(y5(:,:,:,iSet));
            subplot(5,5,iSet);
            plot(pwr1./pwr0,'.-'); hold on;
            if iSet==21; xlabel('component'); ylabel('score'); else; set(gca,'xtick',[]); end
            title(iSet);
            drawnow
        end
    end
end

