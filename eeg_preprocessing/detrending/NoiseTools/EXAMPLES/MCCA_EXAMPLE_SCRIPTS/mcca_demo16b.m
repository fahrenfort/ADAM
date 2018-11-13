% demo code for MCCA

disp('Phase scrambled surrogate data to test significance of MCCA components');

%{
EEG data in response to repeated presentation of a 1 kHz tone for 25 subjects.

Data are organized as a 4D matrix with dimensions time, channel, trial, subject.

Data are preprocessed prior to epoching. Bad channels are detected automatically and removed, 
and data are downsampled to 128 Hz and filtered or detrended (several versions).

There are at least 124 repeats for each subject. We can use this to know
what is the 'true' cortical response to stimulation, and use it a 'ground
truth' to evaluate MCCA.

%}


clear
fnames={
    '../../DATA/MCCA_EXAMPLE_DATA/localiser_nofilt.mat'
	'../../DATA/MCCA_EXAMPLE_DATA/localiser_detrend_10.mat'
	'../../DATA/MCCA_EXAMPLE_DATA/localiser_detrend_epochs_1.mat'
	'../../DATA/MCCA_EXAMPLE_DATA/localiser_detrend_epochs_10.mat'
	'../../DATA/MCCA_EXAMPLE_DATA/localiser_hpf_1Hz.mat'
	'../../DATA/MCCA_EXAMPLE_DATA/localiser_hpf_0.1Hz.mat'
    };
fname=fnames{5};
if exist(fname, 'file')
    load (fname); % loads x sr preStim postStim
else
    error('Download example data from http://audition.ens.fr/adc/NoiseTools/DATA/MCCA_EXAMPLE_DATA/');
end

[nsamples,nchans,ntrials,nsets]=size(x);
t=linspace(-preStim/sr, postStim/sr, nsamples);

%{
Apply DSS to emphasize components repeatable across trials.  The aim is to
reveal the brain activity that we expect MCCA to reveal.
%}


NKEEP=20; % number of PCs to keep in initial dimensionality reduction (optional)
NKEEP2=5; % number of DSS components to keep for each subject

for dimFlag = [0 1] % repeat analysis without and with dimensionality reduction

    figure(1+dimFlag);  clf; %hold on; %clf;
    if dimFlag==0
        nt_banner('no dimensionality reduction');
    else
        nt_banner('dimensionality reduced to 20 per subject');
    end

    % first apply the analysis to the real data (phase scrambled next)
    zz=zeros(nsamples,NKEEP2,ntrials,nsets);
    a=[];
    for iSet=1:nsets
        xx=x(:,:,:,iSet);
        xx=nt_demean(xx,find(t<0)); % baseline correction   
        if 0
            % randomize phase (this should scramble all patterns).
            xx=nt_phase_scramble(xx);
            xx=nt_demean(xx,find(t<0)); % baseline correction
        end 
        if dimFlag
            % reduce dimensionality
            xx=nt_pca(xx);
            xx=xx(:,1:NKEEP,:);
        end
        % apply DSS to emphasize repeatable components
        [todss,pwr0,pwr1]=nt_dss1(xx); 
        a(iSet,:)=pwr1./pwr0;
        z=nt_mmat(xx,todss);    
        % gather selected DSS components of each subjects
        zz(:,:,:,iSet)=z(:,1:NKEEP2,:);
    end
    % apply MCCA to selected DSS components
    zz=permute(zz,[1 3 2 4]); % --> nsamples x ntrials x NKEEP2 x nsets
    zz=reshape(zz,[nsamples*ntrials,NKEEP2,nsets]); % --> nsamples*ntrials x NKEEP2 x nsets
    zz=zz(:,:); % --> nsamples*ntrials x NKEEP2*nsets
    C=zz'*zz;
    [A,score,AA]=nt_mcca(C,NKEEP2);
    zzz=zz*A;
    zzz=reshape(zzz,[nsamples, ntrials, NKEEP2*nsets]);
    zzz=permute(zzz,[1 3 2]);
    plot(score(1:10),'.-', 'linewidth',2); title('SC variance profile'); 
    xlabel('SC'); ylabel('variance');
    hold on


    % repeated trials with phase scramble
    for iRepeat=1:10
        % first apply the analysis to the real data (phase scrambled next)
        zz=zeros(nsamples,NKEEP2,ntrials,nsets);
        for iSet=1:nsets
            xx=x(:,:,:,iSet);
            xx=nt_demean(xx,find(t<0)); % baseline correction   
            if 1
                % randomize phase (this should scramble all patterns).
                xx=nt_phase_scramble(xx);
                xx=nt_demean(xx,find(t<0)); % baseline correction
            end 
            if dimFlag
                % reduce dimensionality
                xx=nt_pca(xx);
                xx=xx(:,1:NKEEP,:);
            end
            % apply DSS to emphasize repeatable components
            [todss,pwr0,pwr1]=nt_dss1(xx); 
            a(iSet,:)=pwr1./pwr0;
            z=nt_mmat(xx,todss);    
            % gather selected DSS components of each subjects
            zz(:,:,:,iSet)=z(:,1:NKEEP2,:);
        end
        % apply MCCA to selected DSS components
        zz=permute(zz,[1 3 2 4]); % --> nsamples x ntrials x NKEEP2 x nsets
        zz=reshape(zz,[nsamples*ntrials,NKEEP2,nsets]); % --> nsamples*ntrials x NKEEP2 x nsets
        zz=zz(:,:); % --> nsamples*ntrials x NKEEP2*nsets
        C=zz'*zz;
        [A,score,AA]=nt_mcca(C,NKEEP2);
        zzz=zz*A;
        zzz=reshape(zzz,[nsamples, ntrials, NKEEP2*nsets]);
        zzz=permute(zzz,[1 3 2]);
        plot(score(1:10),'.-'); title('SC variance profile'); 
        xlabel('SC'); ylabel('variance');
    end
    legend('intact', 'phase scrambled, 10 repeats'); legend boxoff

end

