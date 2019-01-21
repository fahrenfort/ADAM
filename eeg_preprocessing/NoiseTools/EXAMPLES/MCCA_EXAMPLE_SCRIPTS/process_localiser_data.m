% process localiser data

clear
original_sr = 2048; % Hz
sr = 128;
dname='/Users/adc/05/WORK/MCCA/Localiser_data/';
addpath (dname);

if 1
    % read data files, save in matlab format
    for subj = [1:3 5:26]
        disp(subj)

        fileName = [dname, num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.bdf'];
        [EEG_raw, trigs] = Read_bdf(fileName);
        eeg=EEG_raw';

        % Getting triggers
        trigs=trigs-min(trigs);
        trigs(trigs>256) = trigs(trigs>256)-min(trigs(trigs>256));

        % Getting the trials starting points
        startSamples = find(trigs>0);
        startSamples = startSamples(find([1 diff(startSamples)~=1]));
        
        % smooth at 1/50Hz
        eeg=nt_smooth(eeg,original_sr/50);

        % downsample
        eeg = nt_resample(eeg,sr,original_sr);
        startSamples = round(startSamples*sr/original_sr); 

        % eeg = eeg - repmat(mean(eeg,1),128,1);       
        outputFileName = ['/Users/adc/05/WORK/MCCA/Localiser_data/' num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.mat'];
        disp(['Saving ' outputFileName])
        save(outputFileName, 'eeg', 'startSamples', 'sr');
    end
end

if 1
    subjs = [1:3, 5:26];
    for iSubj=1:numel(subjs)
        subj=subjs(iSubj);
        disp(subj);
        fileName = ['/Users/adc/05/WORK/MCCA/Localiser_data/' num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.mat'];
        load(fileName, 'eeg', 'startSamples')
        eeg = eeg(:,1:64);        
        eeg = double(removeBadChannels(eeg')'); % EEGLAB bad-channels rejection
    %     % referencing
    %     eeg = eeg - repmat(mean(eeg,2),1,64);       % referencing
    %     load('LPF30_128Hz.mat')
    %     eeg = filtfilthd(Hd,eeg);
        fs = 128; % Hz
        preStim = round(0.3 * fs); % ms
        postStim = round(0.5 * fs); % ms
        clear epochs
        countEp = 1;
        for idx = startSamples
            epochs(:,:,countEp) = eeg(idx-preStim:idx+postStim,:);
            countEp = countEp + 1;
        end
        epochs=nt_demean(epochs);
        epochsAll{iSubj} = epochs;
    end
    nsubj=numel(epochsAll);
    [nsamples,nchans,ntrials]=size(epochsAll{1});
    for iSubj=1:nsubj
        ntrials=min(ntrials,size(epochsAll{iSubj},3));
    end
    x=zeros(nsamples,nchans,ntrials,nsubj);
    for iSubj=1:nsubj
        x(:,:,:,iSubj)=nt_demean(epochsAll{iSubj}(:,:,1:ntrials));
    end
    save ../../DATA/MCCA_EXAMPLE_DATA/localiser_nofilt x sr preStim postStim
end

if 1
        subjs = [1:3, 5:26];
    for iSubj=1:numel(subjs)
        subj=subjs(iSubj);
        disp(subj);
        fileName = ['/Users/adc/05/WORK/MCCA/Localiser_data/' num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.mat'];
        load(fileName, 'eeg', 'startSamples')
        eeg = eeg(:,1:64);        
        eeg = double(removeBadChannels(eeg')'); % EEGLAB bad-channels rejection
    %     % referencing
    %     eeg = eeg - repmat(mean(eeg,2),1,64);       % referencing
    %     load('LPF30_128Hz.mat')
    %     eeg = filtfilthd(Hd,eeg);
        fs = 128; % Hz
        HPF=1; %Hz
        [B,A]=butter(2,HPF/(fs/2),'high');
        eeg=filter(B,A,eeg);
        preStim = round(0.3 * fs); % ms
        postStim = round(0.5 * fs); % ms
        clear epochs
        countEp = 1;
        for idx = startSamples
            epochs(:,:,countEp) = eeg(idx-preStim:idx+postStim,:);
            countEp = countEp + 1;
        end
        epochs=nt_demean(epochs);
        epochsAll{iSubj} = epochs;
    end
    nsubj=numel(epochsAll);
    [nsamples,nchans,ntrials]=size(epochsAll{1});
    for iSubj=1:nsubj
        ntrials=min(ntrials,size(epochsAll{iSubj},3));
    end
    x=zeros(nsamples,nchans,ntrials,nsubj);
    for iSubj=1:nsubj
        x(:,:,:,iSubj)=nt_demean(epochsAll{iSubj}(:,:,1:ntrials));
    end
    save ../../DATA/MCCA_EXAMPLE_DATA/localiser_hpf_1Hz x sr preStim postStim
end

if 1
        subjs = [1:3, 5:26];
    for iSubj=1:numel(subjs)
        subj=subjs(iSubj);
        disp(subj);
        fileName = ['/Users/adc/05/WORK/MCCA/Localiser_data/' num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.mat'];
        load(fileName, 'eeg', 'startSamples')
        eeg = eeg(:,1:64);        
        eeg = double(removeBadChannels(eeg')'); % EEGLAB bad-channels rejection
    %     % referencing
    %     eeg = eeg - repmat(mean(eeg,2),1,64);       % referencing
    %     load('LPF30_128Hz.mat')
    %     eeg = filtfilthd(Hd,eeg);
        fs = 128; % Hz
        HPF=.1; %Hz
        [B,A]=butter(2,HPF/(fs/2),'high');
        eeg=filter(B,A,eeg);
        preStim = round(0.3 * fs); % ms
        postStim = round(0.5 * fs); % ms
        clear epochs
        countEp = 1;
        for idx = startSamples
            epochs(:,:,countEp) = eeg(idx-preStim:idx+postStim,:);
            countEp = countEp + 1;
        end
        epochs=nt_demean(epochs);
        epochsAll{iSubj} = epochs;
    end
    nsubj=numel(epochsAll);
    [nsamples,nchans,ntrials]=size(epochsAll{1});
    for iSubj=1:nsubj
        ntrials=min(ntrials,size(epochsAll{iSubj},3));
    end
    x=zeros(nsamples,nchans,ntrials,nsubj);
    for iSubj=1:nsubj
        x(:,:,:,iSubj)=nt_demean(epochsAll{iSubj}(:,:,1:ntrials));
    end
    save ../../DATA/MCCA_EXAMPLE_DATA/localiser_hpf_0.1Hz.mat x sr preStim postStim
end

if 1
    subjs = [1:3, 5:26];
    for iSubj=1:numel(subjs)
        subj=subjs(iSubj);
        disp(subj);
        fileName = ['/Users/adc/05/WORK/MCCA/Localiser_data/' num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.mat'];
        load(fileName, 'eeg', 'startSamples')
        eeg = eeg(:,1:64);        
        eeg = double(removeBadChannels(eeg')'); % EEGLAB bad-channels rejection
    %     % referencing
    %     eeg = eeg - repmat(mean(eeg,2),1,64);       % referencing
    %     load('LPF30_128Hz.mat')
    %     eeg = filtfilthd(Hd,eeg);
        fs = 128; % Hz
        ORDER=10;
        eeg=nt_detrend(eeg,ORDER);
        preStim = round(0.3 * fs); % ms
        postStim = round(0.5 * fs); % ms
        clear epochs
        countEp = 1;
        for idx = startSamples
            epochs(:,:,countEp) = eeg(idx-preStim:idx+postStim,:);
            countEp = countEp + 1;
        end
        epochs=nt_demean(epochs);
        epochsAll{iSubj} = epochs;
    end
    nsubj=numel(epochsAll);
    [nsamples,nchans,ntrials]=size(epochsAll{1});
    for iSubj=1:nsubj
        ntrials=min(ntrials,size(epochsAll{iSubj},3));
    end
    x=zeros(nsamples,nchans,ntrials,nsubj);
    for iSubj=1:nsubj
        x(:,:,:,iSubj)=nt_demean(epochsAll{iSubj}(:,:,1:ntrials));
    end
    save ../../DATA/MCCA_EXAMPLE_DATA/localiser_detrend_10 x sr preStim postStim
end

if 1
    subjs = [1:3, 5:26];
    for iSubj=1:numel(subjs)
        subj=subjs(iSubj);
        disp(subj);
        fileName = ['/Users/adc/05/WORK/MCCA/Localiser_data/' num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.mat'];
        load(fileName, 'eeg', 'startSamples')
        eeg = eeg(:,1:64);        
        eeg = double(removeBadChannels(eeg')'); % EEGLAB bad-channels rejection
    %     % referencing
    %     eeg = eeg - repmat(mean(eeg,2),1,64);       % referencing
    %     load('LPF30_128Hz.mat')
    %     eeg = filtfilthd(Hd,eeg);
        fs = 128; % Hz
        EXTRA=1*fs; % extra padding to support detrending
        preStim = round(0.3 * fs)+EXTRA; % ms
        postStim = round(0.5 * fs)+EXTRA; % ms
        clear epochs
        countEp = 1;
        for idx = startSamples
            epochs(:,:,countEp) = eeg(idx-preStim:idx+postStim,:);
            countEp = countEp + 1;
        end
        epochs=nt_demean(epochs);
        ORDER=1;
        epochs=nt_detrend(epochs,ORDER);
        epochs=epochs(EXTRA+1:end-EXTRA,:,:);
        epochsAll{iSubj} = epochs;
    end
    nsubj=numel(epochsAll);
    [nsamples,nchans,ntrials]=size(epochsAll{1});
    for iSubj=1:nsubj
        ntrials=min(ntrials,size(epochsAll{iSubj},3));
    end
    x=zeros(nsamples,nchans,ntrials,nsubj);
    for iSubj=1:nsubj
        x(:,:,:,iSubj)=nt_demean(epochsAll{iSubj}(:,:,1:ntrials));
    end
    save ../../DATA/MCCA_EXAMPLE_DATA/localiser_detrend_epochs_1 x sr preStim postStim
end

if 1
    subjs = [1:3, 5:26];
    for iSubj=1:numel(subjs)
        subj=subjs(iSubj);
        disp(subj);
        fileName = ['/Users/adc/05/WORK/MCCA/Localiser_data/' num2str(subj) '/Localiser_daniel_1_' num2str(subj) '_64.mat'];
        load(fileName, 'eeg', 'startSamples')
        eeg = eeg(:,1:64);        
        eeg = double(removeBadChannels(eeg')'); % EEGLAB bad-channels rejection
    %     % referencing
    %     eeg = eeg - repmat(mean(eeg,2),1,64);       % referencing
    %     load('LPF30_128Hz.mat')
    %     eeg = filtfilthd(Hd,eeg);
        fs = 128; % Hz
        EXTRA=1*fs; % extra padding to support detrending
        preStim = round(0.3 * fs)+EXTRA; % ms
        postStim = round(0.5 * fs)+EXTRA; % ms
        clear epochs
        countEp = 1;
        for idx = startSamples
            epochs(:,:,countEp) = eeg(idx-preStim:idx+postStim,:);
            countEp = countEp + 1;
        end
        epochs=nt_demean(epochs);
        ORDER=10;
        epochs=nt_detrend(epochs,ORDER);
        epochs=epochs(EXTRA+1:end-EXTRA,:,:);
        epochsAll{iSubj} = epochs;
    end
    nsubj=numel(epochsAll);
    [nsamples,nchans,ntrials]=size(epochsAll{1});
    for iSubj=1:nsubj
        ntrials=min(ntrials,size(epochsAll{iSubj},3));
    end
    x=zeros(nsamples,nchans,ntrials,nsubj);
    for iSubj=1:nsubj
        x(:,:,:,iSubj)=nt_demean(epochsAll{iSubj}(:,:,1:ntrials));
    end
    save ../../DATA/MCCA_EXAMPLE_DATA/localiser_detrend_epochs_10 x sr preStim postStim
end
