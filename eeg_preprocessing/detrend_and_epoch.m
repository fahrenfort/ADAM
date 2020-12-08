function detrend_and_epoch(datadir,filename,outputdir, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, mask_bad_data, channelpool, remove_bad_chans, event_codes)
% function detrend_and_epoch(datadir,filename,outputdir, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, mask_bad_data, channelpool, remove_bad_chans, event_codes)
% detrend_and_epoch is an internal function of the ADAM toolbox. Loads EEGLAB data and performs a
% multivariate classification procedure. Refer to the help of ADAM_DETREND_AND_EPOCH for proper
% instructions on how to use this function.
%
% Internal function of the ADAM toolbox by J.J.Fahrenfort, VU 2014, 2015, 2018, 2019, 2020
%
% See also: ADAM_DETREND_AND_EPOCH

%% input checks, convert inputs to doubles etc 
if nargin < 12
    help detrend_and_epoch;
    error('not enough input arguments');
end
if ischar(start_epoch)
    start_epoch = string2double(start_epoch);
end
if ischar(end_epoch)
    end_epoch = string2double(end_epoch);
end
if ischar(polynomial_order)
    polynomial_order = string2double(polynomial_order);
end
if ischar(pad_length)
    pad_length = string2double(pad_length);
end
if ischar(start_mask)
    start_mask = string2double(start_mask);
end
if ischar(end_mask)
    end_mask = string2double(end_mask);
end
if ischar(mask_only_current)
    if strcmpi(mask_only_current,'yes')
        preset_mask_on_trial = 'current';
    else
        preset_mask_on_trial = 'all';
    end
    if start_mask >= end_mask
        preset_mask_on_trial = 'none';
    end
end
if ischar(mask_bad_data)
    if strcmpi(mask_bad_data,'yes')
        mask_bad_data = true;
    else
        mask_bad_data = false;
    end
end
if ischar(remove_bad_chans)
    if strcmpi(remove_bad_chans,'yes')
        remove_bad_chans = true;
    else
        remove_bad_chans = false;
    end
end
% some checks
if isempty(end_mask)
    end_mask = end_epoch - 0.25;
    disp('warning: no end point for mask was given, assuming 250 ms before end of trial for end mask');
end
if isempty(start_mask)
    start_mask = 0;
    disp('warning: no starting point for mask was given, assuming time = 0 for start mask');
end
if isempty(pad_length)
    pad_length = 50;
end
if isempty(polynomial_order)
    polynomial_order = 30;
end
if isempty(outputdir)
    outputdir = datadir;
else
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
end
% FIX TIME, INTERNALLY EVERYTHING IS IN MILLISECONDS
start_epoch = start_epoch * 1000;
end_epoch = end_epoch * 1000;
start_mask = start_mask * 1000;
end_mask = end_mask * 1000;
% conditions should be defined as a string of comma separated values
conditions = string2double(event_codes);

[~,fname,~] = fileparts(filename);

%% load EEGLAB data and make sure there is electrode position data
EEG = pop_loadset('filename',[fname '.set'],'filepath',datadir);
% next identify bad channels
try
    eeg_channels = select_channels({EEG.chanlocs(:).labels},channelpool);
catch
    error('Stopping now, there are no EEG channels in this set? ADAM only works with standard 10-20 EEG labels.');
end
% double check whether channel location information is present
nopos_channels = [];
for cEl=1:length(EEG.chanlocs)
    if (any(isempty(EEG.chanlocs(1,cEl).X)&isempty(EEG.chanlocs(1,cEl).Y)&isempty(EEG.chanlocs(1,cEl).Z)&isempty(EEG.chanlocs(1,cEl).theta)&isempty(EEG.chanlocs(1,cEl).radius)))
        nopos_channels = [nopos_channels cEl];
    end
end
EEG = pop_select(EEG, 'channel', eeg_channels);
% look up electrode info
if any(ismember(eeg_channels,nopos_channels))
    disp(['WARNING: Channels ' num2str(nopos_channels) ' have incomplete location information. Now attempting to read in location information.']);
    EEG = pop_chanedit(EEG, 'lookup', trycapfile);
end

%% identify bad channels and bad data
if remove_bad_chans || mask_bad_data
    % strong 1 Hz highpass filter (temporary, only for identifying bad channels/data)
    EEG_filt = pop_eegfiltnew(EEG, 1);
end
if remove_bad_chans
    % identify and remove bad channels
    EEG_nobadchans = clean_channels(EEG_filt);
    orig_chanlocs = EEG.chanlocs;
    clean_chanlocs = EEG_nobadchans.chanlocs;
    rej_channels = setdiff({orig_chanlocs.labels},{clean_chanlocs.labels});
    % reject bad channels
    EEG = pop_select(EEG, 'nochannel', rej_channels);
    EEG_filt = pop_select(EEG_filt, 'nochannel', rej_channels);
end
if mask_bad_data 
    % EEG_filt = pop_eegfiltnew(EEG_filt, 110, 140); % interested in muscle artefacts in particular, so band-pass between 110 and 140
    % create mask of the data that contains clean data (to mask out dirty parts for detrending)
    for cChan = 1:size(EEG_filt.data,1)
        temp_EEG = pop_select(EEG_filt, 'channel', cChan);
        [~,clean_mask(cChan,:)] = clean_windows(temp_EEG,[],[-25,25]);
    end
    %[~,clean_mask] = clean_windows(EEG_filt);
end
clear EEG_filt EEG_nobadchans temp_EEG;

%% filter, detrend or nothing

% little hack: apply filtering if polynomial_order < 0, so we can compare filtering to detrending
if polynomial_order < 0
    cutoff = - polynomial_order;
    ncycles = 3;
    filter_order = ncycles.*(EEG.srate/cutoff);
    filter_order = ceil(filter_order ./ 2) .* 2;
    % windowed sinc FIR filter with Kaiser window type (Widmann et al 2015 JNM)
    EEG = pop_firws(EEG,'fcutoff',cutoff,'forder',filter_order,'ftype','highpass','wtype','kaiser','warg',pop_kaiserbeta(.001));
end

% convert to FT_EEG format for ease of processing
FT_EEG = eeglab2ft(EEG,[],true,[start_epoch end_epoch]/1000); % third argument indicated data is continuous, fourth argument is to indicate epoch window in seconds
clear EEG;

% obtain names to work with
eeg_data = FT_EEG.trial{1};
% FT_EEG is in seconds, but we want eeg_time in milliseconds, so converted here
eeg_time = FT_EEG.time{1}*1000; % FIX TIME, INTERNALLY EVERYTHING IS IN MILLISECONDS
srate = FT_EEG.fsample;
trialinfo = FT_EEG.trialinfo;
label = FT_EEG.label;
trialinfo(:,2) = trialinfo(:,2); % in milliseconds

% remove mean from every channel % WHY???? DON'T DO THAT!!!
% eeg_data = eeg_data-repmat(mean(eeg_data,2),[1 size(eeg_data,2)]);

% mirror-pad edges of the unepoched data, so that extracting wide padded epochs will not be problematic
eeg_data = padarray(eeg_data,[0 pad_length*srate],'both','symmetric');

% eeg_time_old = padarray(eeg_time,[0 pad_length*srate],NaN,'both'); -> this fills it with NaNs which we don't want
% instead, let's mirror-pad the time array with time (also going negative) 
eeg_time_step = (eeg_time(end)-eeg_time(1))/(numel(eeg_time)-1); % determine step size
eeg_time = (eeg_time(1)-(eeg_time_step*pad_length*srate)):eeg_time_step:(eeg_time(end)+(eeg_time_step*pad_length*srate)); % create time line

% create a mask for all trials
if mask_bad_data
    eeg_mask = padarray(clean_mask,[0 pad_length*srate],'both','symmetric'); % pad the clean_mask
else
    eeg_mask = ones(size(eeg_data)); % no masking of bad data, just use the mirror-padded eeg_data
end

% select relevant events
trialinfo = trialinfo(ismember(trialinfo(:,1),conditions),:);

% mask out all trials: do this when wanting to mask out all trials during detrending (instead of only the current trial) 
if strcmpi(preset_mask_on_trial, 'all')
    for cTrials = 1:size(trialinfo,1)
        mask_startind = nearest(eeg_time,trialinfo(cTrials,2)+start_mask);
        mask_stopind = nearest(eeg_time,trialinfo(cTrials,2)+end_mask);
        eeg_mask(:,mask_startind:mask_stopind) = 0; % bugfix, missing channel in mask
    end
end

% initialize trial array
new_trial = NaN([size(trialinfo,1) numel(label) numel(nearest(eeg_time,0):nearest(eeg_time,(end_epoch-start_epoch)))]);
if polynomial_order > 0
    old_trial = new_trial;
end

%% Do the epoching 
for cTrials = 1:size(trialinfo,1)
    
    % identify zero point and borders of trial
    start_ind = nearest(eeg_time,trialinfo(cTrials,2)+start_epoch);
    zero_ind = nearest(eeg_time,trialinfo(cTrials,2));
    stop_ind = nearest(eeg_time,trialinfo(cTrials,2)+end_epoch);
    
    % extract wide padded trial
    pad_time = eeg_time((start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
    pad_data = eeg_data(:,(start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
    
    temp_mask = eeg_mask; % copy from the source
    if strcmpi(preset_mask_on_trial, 'current') % mask only the current trial
        mask_startind = nearest(eeg_time,trialinfo(cTrials,2)+start_mask);
        mask_stopind = nearest(eeg_time,trialinfo(cTrials,2)+end_mask);
        temp_mask(:,mask_startind:mask_stopind) = 0; % bugfix, missing channel in mask
    end
    % when preset_mask_on_trial == 'none', no preset masking is be applied
    pad_mask = temp_mask(:,(start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate)); % bugfix, missing channel in mask
    clear temp_mask;
    
    % Do detrending, ONLY if polynomial_order > 0
    if polynomial_order > 0
        % estimate and subtract polynomial
        disp(['Polynomial detrending trial ' num2str(cTrials) ' of ' num2str(size(trialinfo,1)) '...']);
        % create a mask matrix
        wt = pad_mask; % bugfix, already had channel in mask
        % make sure data is physically removed before calling nt_detrending function, just to be sure
        pad_data_orig = pad_data;
        % no way in hell that the detrending function can make use of brain data:
        pad_data_nocogdata = pad_data .* wt;
        % call nt_detrend without brain data
        [tmp,w1,~,regressline1] = nt_detrend(pad_data_nocogdata',1,wt'); % start with 1st order
        % take out the regression line manually
        pad_data = pad_data_orig - regressline1';
        % call nt_detrend with higher order polynomial
        [~, w2,~,regressline2] = nt_detrend(tmp,polynomial_order,w1); % then nth order with mask of previous step
        % take out the regression line manually
        clean_data = pad_data - regressline2'; % manually take out regression slope from actual data
        % epoch original data to a narrow window
        old_trial(cTrials,:,:) = pad_data_orig(:,(pad_length/2*srate+1):end-pad_length/2*srate);
    else
        clean_data = pad_data;
    end
    
    % epoch data to a narrow window
    new_trial(cTrials,:,:) = clean_data(:,(pad_length/2*srate+1):end-pad_length/2*srate);
    new_time = pad_time - eeg_time(zero_ind);
    trial_time = new_time((pad_length/2*srate+1):end-pad_length/2*srate);
    
    %% plot a trial halfway through the experiment for illustration purposes
    if cTrials == round(size(trialinfo,1)/2) && polynomial_order > 0
    % if cTrials == 8 && polynomial_order > 0 % for figure 3, applied to subject 2, remove
        % create an invisible image, to save later on
        f = figure('units','normalized','outerposition',[0 0 1 1],'visible', 'off');
        % f = figure('units','normalized','outerposition',[0 0 1 1],'visible', 'on'); % for figure 3, remove
        % some settings
        sample_index = round(linspace(1,size(pad_data_orig,2),1000)); % (pad_length/2*srate+1):end-pad_length/2*srate
        trial_index = find(new_time > start_epoch & new_time < end_epoch);
        EEG_index = select_channels(label,'EEG'); % identify the electrodes with the largest standard deviation
        [~, order] = sort(std(pad_data_orig'),'descend');
        order = order(ismember(order,EEG_index));
        % make mask for first 5 electrodes, top plot
        elec_index = order(1:5);
        % elec_index = order(1:3); % for figure 3, remove
        plot_mask = w1;
        plot_mask(plot_mask==0) = NaN;
        plot_mask = plot_mask(:,elec_index);
        top = max(max(pad_data_orig(elec_index,sample_index)));
        for c = 1:numel(elec_index)
            plot_mask(plot_mask(:,c)==1,c) = top + (top/5) - c*(top/75);
        end
        % top plot
        subplot(4,2,[1 2]);
        hold on;
        plot(new_time(sample_index)/1000,pad_data_orig(elec_index,sample_index)'); % original data
        plot(new_time(sample_index)/1000,regressline1(sample_index,elec_index(1)),'k--');
        plot(new_time(sample_index)/1000,plot_mask(sample_index,:),'k','LineWidth',1);
        plot(new_time(sample_index)/1000,regressline1(sample_index,elec_index(2:end)),'k--');
        legend({label{elec_index} '1st order' 'masks'}); %' ' ' ' ['mask ' labels{1}] 'mask c3' 'mask f1'
        title('raw data');
        minlim = min(min(pad_data_orig(elec_index,sample_index)));
        maxlim = top+(top/5);
        ylim([minlim-top/5 maxlim]);
        % ylim([-3 3]); % for figure 3, remove
        xlim([min(new_time) max(new_time)+(max(new_time)-min(new_time))/4]/1000);
        xlim([-25 31.5]); % for figure 3, remove
        % make mask for 10 electrodes, middle plot
        plot_mask = w2;
        plot_mask(plot_mask==0) = NaN;
        plot_mask = plot_mask(:,elec_index);
        top = max(max(pad_data(elec_index,trial_index)));
        for c = 1:numel(elec_index)
            plot_mask(plot_mask(:,c)==1,c) = top + (top/5) - c*(top/75);
        end
        % middle plot
        subplot(4,2,[3 4]);
        hold on;
        plot(new_time(sample_index)/1000,pad_data(elec_index,sample_index)'); % first polynomial removed
        plot(new_time(sample_index)/1000,regressline2(sample_index,elec_index(1)),'k--');
        plot(new_time(sample_index)/1000,plot_mask(sample_index,:),'k','LineWidth',1);
        plot(new_time(sample_index)/1000,regressline2(sample_index,elec_index(2:end)),'k--');
        legend({label{elec_index} [num2str(polynomial_order) 'th order'] 'masks'});
        title('after 1st order polynomial removal');
        minlim = min(min(pad_data(elec_index,trial_index)));
        maxlim = top+(top/5);
        ylim([minlim+minlim/10 maxlim]);
        % ylim([-2.5 2.5]); % for figure 3, remove
        xlim([min(new_time) max(new_time)+(max(new_time)-min(new_time))/4]/1000);
        % xlim([-25 31.5]); % for figure 3, remove
        % bottom plot
        subplot(4,2,[5 6]);
        plot(new_time(sample_index)/1000,clean_data(elec_index,sample_index)'); % final
        xlabel('time in seconds');
        legend({label{elec_index}});
        title(['after ' num2str(polynomial_order) 'th order polynomial removal']);
        minlim = min(min(clean_data(elec_index,trial_index)));
        maxlim = max(max(clean_data(elec_index,trial_index)));
        ylim([minlim+minlim/10  maxlim+maxlim/10]);
        % ylim([-1 1]); % for figure 3, remove
        xlim([min(new_time) max(new_time)+(max(new_time)-min(new_time))/4]/1000);
        % xlim([-25 31.5]); % for figure 3, remove
        % print -painters -dsvg '/Users/VU-MBP/Dropbox/Work/Artikelen/- Collegas - coauteur/Joram van Driel/Filtering_paper/JoN methods/figures/orig/figure3.svg'; % remove
    end
end

%% Plot butterfly plots of ERP. Only plot EEG channels (no ocular channels)
% Also baseline < 0 (baseline only for illustration purposes, not in saved data)
if polynomial_order > 0
    oldERP = mean(old_trial(:,EEG_index,:) - repmat(mean(old_trial(:,EEG_index,trial_time<0),3),[1 1 size(old_trial,3)]),1);
    clear old_trial;
    newERP = mean(new_trial(:,EEG_index,:) - repmat(mean(new_trial(:,EEG_index,trial_time<0),3),[1 1 size(new_trial,3)]),1);
    subplot(4,2,7);
    plot(trial_time, squeeze(oldERP));
    xlabel('time in milliseconds');
    ntitle('ERP butterfly before polynomial removal','FontSize',8,'FontWeight','bold');
    xlim([min(trial_time) max(trial_time)]);
    subplot(4,2,8);
    plot(trial_time, squeeze(newERP));
    xlabel('time in milliseconds');
    ntitle('ERP butterfly after polynomial removal','FontSize',8,'FontWeight','bold');
    xlim([min(trial_time) max(trial_time)]);
    print(f,'-dpng', [outputdir filesep fname '.png']);
    close(f);
end

%% create a new epoched and detrended FT_EEG dataset and convert back to EEGLAB format
DETRENDED_FT_EEG.dimord = 'rpt_chan_time';
DETRENDED_FT_EEG.label = label;
DETRENDED_FT_EEG.time = trial_time; % in milliseconds
DETRENDED_FT_EEG.trial = new_trial; 
DETRENDED_FT_EEG.trialinfo = trialinfo(:,1);
DETRENDED_FT_EEG.fsample = srate;

EEG = ft2eeglab(DETRENDED_FT_EEG);

%% now that detrending is complete, interpolate faulty electrodes if required
if remove_bad_chans && ~isempty(rej_channels)
    EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
    fid = fopen([outputdir filesep 'bad_channels_' fname '.txt'], 'wt' );
    for c = 1:numel(rej_channels)
        fprintf( fid, '%s\n', rej_channels{c});
    end
    fclose(fid);
end

%% save
% can also save as FT_EEG if that is preferred:
% save([outputdir filesep fname '.mat'],'DETRENDED_FT_EEG'); 
% OR as EEGLAB file
pop_saveset(EEG, 'filename',[fname '.set'],'filepath',outputdir);