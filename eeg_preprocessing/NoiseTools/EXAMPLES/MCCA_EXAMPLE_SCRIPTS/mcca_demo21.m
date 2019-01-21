% demo code for MCCA

disp('EEG data: the importance of preprocessing');

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
banners={
    'no filtering'
    'detrend raw with polynomial order 10'
    'detrend each epoch with polynomial order 1'
    'detrend each epoch with polynomial order 10'
    'highpass raw butter order 2, HPF=1 Hz'
    'highpass raw butter order 2, HPF=0.1 Hz'
    };

for iFname=[1 4 5 6] %1:numel(fnames);
    fname=fnames{iFname};
    if exist(fname, 'file')
        load (fname); % loads x sr preStim postStim 
    else
        error('Download example data from http://audition.ens.fr/adc/NoiseTools/DATA/MCCA_EXAMPLE_DATA/');
    end

    [nsamples,nchans,ntrials,nsets]=size(x);
    t=linspace(-preStim/sr, postStim/sr, nsamples);

    % remove mean over pre-stimulus segment ('baseline correction')
    for iSet=1:nsets
        x(:,:,:,iSet)=nt_demean(x(:,:,:,iSet),find(t<0));
    end

    % reduce dimensionality
    NKEEP=20; 
    y=zeros(nsamples,NKEEP,ntrials,nsets);
    for iSet=1:nsets
        a=nt_pca(x(:,:,:,iSet)); 
        a=a(:,1:NKEEP,:);
        y(:,:,:,iSet)=a;
    end

    % unfold each data matrix (concatenate trials)
    yy=permute(y,[1 3 2 4]);
    yy=reshape(yy,[nsamples*ntrials,NKEEP,nsets]);

    % concatenate over channels
    yy=yy(:,:);

    % MCCA
    C=yy'*yy;
    [A,score,AA]=nt_mcca(C,NKEEP);
    z=yy*A;

    % fold SCs (--> samples * SCs * trials)
    z=nt_fold(z,nsamples);

    figure(iFname);
    nt_banner(banners{iFname});
    
    subplot 221;
    plot(score, '.-'); title('variance profile'); xlabel('SC'); ylabel('variance');
    subplot 223
    nt_imagescc(nt_normcol(mean(z(:,1:30,:),3))');
    title('first 30 SCs, trial avg');
    ylabel('SC'); xlabel('time (samples)');
    for k=1:3
        subplot (3,2,2*k);
        nt_bsplot(z(:,k,:),[],[],t);
        title(k); 
        set(gca,'ytick',[]); 
    end
    xlabel('time re stim onset (s)');
    drawnow
end
