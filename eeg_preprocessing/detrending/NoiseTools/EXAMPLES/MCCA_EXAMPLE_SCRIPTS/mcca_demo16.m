% demo code for MCCA

disp('Real EEG data: apply DSS to emphasize repeatable activity');

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
    whos
else
    error('Download example data from http://audition.ens.fr/adc/NoiseTools/DATA/MCCA_EXAMPLE_DATA/');
end

[nsamples,nchans,ntrials,nsets]=size(x);
t=linspace(-preStim/sr, postStim/sr, nsamples);

%{
Apply DSS to emphasize components repeatable across trials.  The aim is to
reveal the brain activity that we expect MCCA to reveal.
%}

NKEEP2=5;
zz=zeros(nsamples,NKEEP2,ntrials,nsets);
figure(1); clf; nt_banner('first DSS component, all subjects')
figure(2); clf; nt_banner('first 6 DSS components, all subjects');
for iSet=1:nsets
    xx=x(:,:,:,iSet);
    xx=nt_demean(xx,find(t<0)); % baseline correction
    
    if 0
        % identify, remove outlier trials
        idx=nt_find_outlier_trials(nt_demean2(xx),3);
        xx=nt_demean(xx(:,:,idx),find(t<0));
        idx=nt_find_outlier_trials(nt_demean2(xx),2);
        xx=nt_demean(xx(:,:,idx),find(t<0));
    end
    
    if 0
        % reduce dimensionality
        xx=nt_pca(xx);
        NKEEP=20;
        xx=xx(:,1:NKEEP,:);
    end
    
    if 0
        % randomize phase (this should scramble all patterns).
        xx=nt_phase_scramble(xx);
        xx=nt_demean(xx,find(t<0)); % baseline correction
    end
 
    % apply DSS to emphasize repeatable components
    [todss,pwr0,pwr1]=nt_dss1(xx); 
    a(iSet,:)=pwr1./pwr0;
    z=nt_mmat(xx,todss);
    
    figure(1);
    subplot(5,5,iSet);
    nt_bsplot(z(:,1,:)*sign(mean(mean(z(:,1,:),3))),[],[],t); 
    set(gca,'ytick',[]); title(iSet);
    if iSet ~= 21; set(gca,'xtick', [], 'ytick', []); else; xlabel('time (s)'); end
    
    figure(2);
    subplot(5,5,iSet);
    b=mean(z(:,1:6,:),3);
    b=bsxfun(@times,b,sign(mean(b))); % flip signs to look good
    nt_imagescc(b'); 
    title(iSet); 
    if iSet ~= 21; set(gca,'xtick', [], 'ytick', []); else; xlabel('sample'); ylabel('component'); end
    
    % gather selected DSS components of each subjects
    zz(:,:,:,iSet)=z(:,1:NKEEP2,:);
end

figure(3); clf;
plot(a(:,1:20)', '.-'); ylabel('score'); xlabel('component'); title('DSS analysis for each subject');

% apply MCCA to selected DSS components

zz=permute(zz,[1 3 2 4]); % --> nsamples x ntrials x NKEEP2 x nsets
zz=reshape(zz,[nsamples*ntrials,NKEEP2,nsets]); % --> nsamples*ntrials x NKEEP2 x nsets
zz=zz(:,:); % --> nsamples*ntrials x NKEEP2*nsets
C=zz'*zz;
[A,score,AA]=nt_mcca(C,NKEEP2);
zzz=zz*A;
zzz=reshape(zzz,[nsamples, ntrials, NKEEP2*nsets]);
zzz=permute(zzz,[1 3 2]);

figure(4);  clf; %hold on; %clf;
plot(score,'.-'); title('SC variance profile'); 
xlabel('SC'); ylabel('variance');

figure(5); clf
nt_banner('MCCA applied to selected DSS components');
for k=1:6
    subplot(2,3,k);
    nt_bsplot(zzz(:,k,:),[],[],t);
    title(['SC ',num2str(k)]); set(gca,'ytick',[]);
end
xlabel('time (s)');

