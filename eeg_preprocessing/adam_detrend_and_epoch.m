function adam_detrend_and_epoch(filepath,filenames,outpath, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, varargin)
% function detrend_eeg_and_epoch(filepath,filenames,outpath, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, varargin)
% Reads in CONTINUOUS (!) EEGLAB data, removes faulty electrodes, detrends using polynomials, epochs
% and writes out as EPOCHED (!) EEGLAB file. Note that the data are high-pass filtered early on to
% identify faulty electrodes, but that the highpass filtered data is discarded after identifying
% faulty electrodes. Polynomial detrending is applied to remove drifts in raw data. Faulty
% electrodes are interpolated back in after detrending. filenames: names of files either in a cell
% array or as comma separated string. Wildcards * and ? can be used, e.g. filenames = '*' will take
% all the .set files in the input filepath as sources. The output of this function can serve as
% input for ADAM_MVPA_FIRSTLEVEL.
%
%       start_epoch             in seconds (relative to target events)
%       end_epoch               in seconds (relative to target events) 
%       polynomial_order        order of the polynomial that is used for detrending (Default: 30,
%                               but note that higher orders can improve the data even more).
%       pad_length              wide epoch padding window used to fit the polynomial on, e.g.
%                               pad_length of 100 creates 50 sec pads around both sides of a trial
%                               (Default: 50). This padding window is only used during detrending
%                               and is discarded after the trial itsel is epoched out (when doing
%                               the final epoching).
%       start_mask              in seconds, specifies when your 'cognitive' or other events of
%                               interest start
%       end_mask                in seconds, specifies when your 'cognitive' or other events of
%                               interest end 
%       mask_only_current       false or true (default: true). Specifies whether the algorithm for
%                               any given trial masks out all other trials from the wide window pad,
%                               or only the currently relevant trial. If you have sufficiently long
%                               intertrial intervals you may also set this to false.
%       varargin                a list of conditions to epoch, either as strings or numbers
%                               separated by commas, or as a cell array Use  cond_string([1001:1024
%                               1101:1124]) when called from create_qsub_files in order to epoch all
%                               conditions running from 1001 to 1024 and from 1101 to 1124.
%
% The function also saves a plot for each subject, containg (from top to bottom):
%       1.      A graph of the middle trial in the experiment, plotting the five electrodes with the
%               largest standard deviation across the wide epoch, their 1st order polynomial fits as
%               dotted lines, plus at the top the masks for these five electrodes
%       2.      The same five electrodes, but now with their 1st order polynomial fits removed,
%               their nth order polynomials (n is specified by polynomial_order) as dotted lines,
%               plus at the top the masks for these five electrodes
%       3.      The same five electrodes with their nth order polynomials removed
%       4.      The butterfly ERP before polynomial removal (baselined at <0 for illustration
%               purposes)
%       5.      The butterfly ERP after polynomial removal (baselined at <0 for illustration
%               purposes, the output data of the function are not baselined)
%
% example: adam_detrend_eeg_and_epoch('c:\alldata\','subject01', 'c:\filtereddata', -.5, 1, 30, 100, 0, .6, true, cond_string([1001:1024 1101:1124]));
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2018, 2019
%
% See also ADAM_MVPA_FIRSTLEVEL, ADAM_COMPUTE_GROUP_MVPA, ADAM_PLOT_MVPA

% some input checking
if nargin < 10 || isempty(mask_only_current)
    mask_only_current = true;
end
if nargin < 9 || isempty(end_mask)
    end_mask = end_epoch - 0.25;
    disp('warning: no end point for mask was given, assuming 250 ms before end of trial for end mask');
end
if nargin < 8 || isempty(start_mask)
    start_mask = 0;
    disp('warning: no starting point for mask was given, assuming time = 0 for start mask');
end
if nargin < 7 || isempty(pad_length)
    pad_length = 50;
end
if nargin < 6 || isempty(polynomial_order)
    polynomial_order = 30;
end
if nargin < 5
    error('not enough arguments, need at least filepath, filenames, outpath, start_epoch and end_epoch');
end
if isempty(outpath)
    outpath = filepath;
else
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
end
% convert inputs to doubles
if ischar(start_epoch);
    start_epoch = string2double(start_epoch);
end
if ischar(end_epoch)
    end_epoch = string2double(end_epoch);
end
if ischar(polynomial_order);
    polynomial_order = string2double(polynomial_order);
end
if ischar(pad_length)
    pad_length = string2double(pad_length);
end
if ischar(start_mask);
    start_mask = string2double(start_mask);
end
if ischar(end_mask)
    end_mask = string2double(end_mask);
end
if ischar(mask_only_current)
    if strcmpi(mask_only_current,'yes')
        mask_only_current = true;
    else
        mask_only_current = false;
    end
end
% fix time, internally everything is in milliseconds
start_epoch = start_epoch * 1000;
end_epoch = end_epoch * 1000;
start_mask = start_mask * 1000;
end_mask = end_mask * 1000;
% fix filename
if ~iscell(filenames) && (~isempty(strfind(filenames,'*')) || ~isempty(strfind(filenames,'?')))
    if ~strcmp(filenames(end-3:end),'.set')
        filenames = [filenames '.set'];
    end
    filenames = dir([filepath filesep filenames]);
    filenames = {filenames(:).name};
end
if ~iscell(filenames)
    filenames = regexp(filenames, ',', 'split');
end
% conditions can either be defined as a string of comma separated values,
% or a cell array of numbers or strings, and is converted to a ell array of
% strings to be fed into pop_epoch
if numel(varargin) == 1
    if ischar(varargin{1})
        conditions = regexp(varargin{1}, ',', 'split');
    end
else
    if iscell(varargin)
        conditions = regexp(vec2str([varargin{:}]), ',', 'split');
    end
end
if ~isnumeric(conditions)
    conditions = string2double(conditions);
end

% go
for filename = filenames
    [~,fname,~] = fileparts(filename{1});
    
    %% load EEGLAB data and temporarily apply liberal 1 Hz highpass filter to identify bad channels
    EEG = pop_loadset('filename',[fname '.set'],'filepath',filepath);
    % strong 1 Hz highpass filter
    EEG = pop_eegfiltnew(EEG, 1);
    % find EEG channels, electrode rejection only on EEG
    try 
        eeg_channels = select_channels({EEG.chanlocs(:).labels},'EEG');
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
    % identify and remove bad channels from the EEG, interpolate and write bad channels to text file
    orig_chanlocs = EEG.chanlocs;
    EEG = clean_artifacts(EEG,'Highpass','off','WindowCriterion','off','burst_crit','off');
    clean_chanlocs = EEG.chanlocs;
    rejected_electrodes = setdiff({orig_chanlocs.labels},{clean_chanlocs.labels});
    
    %% load original EEGLAB data back in
    EEG = pop_loadset('filename',[fname '.set'],'filepath',filepath);
    
    % reject bad channels
    EEG = pop_select(EEG, 'nochannel', rejected_electrodes);
    
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
    FT_EEG = eeglab2ft(EEG,[],true); % third argument indicated data is continuous
    clear EEG;
    
    % obtain names to work with
    eeg_data = FT_EEG.trial{1};   
    eeg_time = FT_EEG.time{1}; % in milliseconds
    srate = FT_EEG.fsample;
    trialinfo = FT_EEG.trialinfo; 
    label = FT_EEG.label;
    trialinfo(:,2) = trialinfo(:,2); % in milliseconds
    
    % remove mean from every channel
    eeg_data = eeg_data-repmat(mean(eeg_data,2),[1 size(eeg_data,2)]);
    
    % select relevant events
    trialinfo = trialinfo(ismember(trialinfo(:,1),conditions),:);
    
    % mirror-pad edges of the unepoched data, so that extracting wide padded epochs will not be problematic
    eeg_data = padarray(eeg_data,[0 pad_length*srate],'both','symmetric');
    eeg_time = padarray(eeg_time,[0 pad_length*srate],NaN,'both');
        
    % create a mask for all trials
    eeg_mask = ones(size(eeg_data));
    if ~mask_only_current
        for cTrials = 1:size(trialinfo,1)
            mask_startind = nearest(eeg_time,trialinfo(cTrials,2)+start_mask);
            mask_stopind = nearest(eeg_time,trialinfo(cTrials,2)+end_mask);
            eeg_mask(mask_startind:mask_stopind) = 0;
        end
    end
    
    % do the epoching (and detrending if indicated)
    for cTrials = 1:size(trialinfo,1)
                
        % identify zero point and borders of trial
        start_ind = nearest(eeg_time,trialinfo(cTrials,2)+start_epoch);
        zero_ind = nearest(eeg_time,trialinfo(cTrials,2));
        stop_ind = nearest(eeg_time,trialinfo(cTrials,2)+end_epoch);
        
        % extract wide padded trial
        pad_time = eeg_time((start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
        pad_data = eeg_data(:,(start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
        if mask_only_current % mask only current trial
            pad_mask = ones(1,size(pad_data,2));
            mask_zero_ind = nearest(pad_time,trialinfo(cTrials,2)+start_mask);
            mask_stop_ind = nearest(pad_time,trialinfo(cTrials,2)+end_mask);
            pad_mask(mask_zero_ind:mask_stop_ind) = 0;
        else % mask all trials
            pad_mask = eeg_mask((start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
        end
        
        if polynomial_order > 0
            % estimate and subtract polynomial
            disp(['Polynomial detrending trial ' num2str(cTrials) '...']);
            % create a mask matrix
            wt = repmat(pad_mask,[numel(label) 1]);
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
            clean_data = pad_data-regressline2'; % manually take out regression slope from actual data
            % epoch original data to a narrow window
            old_trial(cTrials,:,:) = pad_data_orig(:,(pad_length/2*srate+1):end-pad_length/2*srate);            
        else
            clean_data = pad_data;
        end
        
        %% epoch data to a narrow window
        new_trial(cTrials,:,:) = clean_data(:,(pad_length/2*srate+1):end-pad_length/2*srate);
        new_time = pad_time - eeg_time(zero_ind);
        trial_time = new_time((pad_length/2*srate+1):end-pad_length/2*srate);
        
        %% plot a trial halfway through the experiment for illustration purposes
        if cTrials == round(size(trialinfo,1)/2) && polynomial_order > 0
            % create an invisible image, to save later on
            f = figure('units','normalized','outerposition',[0 0 1 1],'visible', 'off');
            % some settings
            sample_index = round(linspace(1,size(pad_data_orig,2),1000)); % (pad_length/2*srate+1):end-pad_length/2*srate
            trial_index = find(new_time > start_epoch & new_time < end_epoch);
            EEG_index = select_channels(label,'EEG'); % identify the 10 electrodes with the largest standard deviation
            [~, order] = sort(std(pad_data_orig'),'descend');
            order = order(ismember(order,EEG_index));
            % make mask for first 5 electrodes, top plot
            elec_index = order(1:5);
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
            xlim([min(new_time) max(new_time)+(max(new_time)-min(new_time))/4]/1000);
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
            xlim([min(new_time) max(new_time)+(max(new_time)-min(new_time))/4]/1000);
            % bottom plot
            subplot(4,2,[5 6]);
            plot(new_time(sample_index)/1000,clean_data(elec_index,sample_index)'); % final
            xlabel('time in seconds');
            legend({label{elec_index}});
            title(['after ' num2str(polynomial_order) 'th order polynomial removal']);
            minlim = min(min(clean_data(elec_index,trial_index)));
            maxlim = max(max(clean_data(elec_index,trial_index)));
            ylim([minlim+minlim/10  maxlim+maxlim/10]);
            xlim([min(new_time) max(new_time)+(max(new_time)-min(new_time))/4]/1000);
        end
    end
    
    % Plot butterfly plots of ERP 
    % Only plot EEG channels (no ocular channels), baseline < 0 (baseline only for illustration purposes)
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
        print(f,'-dpng', [outpath filesep fname '.png']);
        close(f);
    end
    
    % create a new epoched and detrended FT_EEG dataset
    DETRENDED_FT_EEG.dimord = 'rpt_chan_time';
    DETRENDED_FT_EEG.label = label;
    DETRENDED_FT_EEG.time = trial_time; % in milliseconds
    DETRENDED_FT_EEG.trial = new_trial;
    DETRENDED_FT_EEG.trialinfo = trialinfo(:,1);
    DETRENDED_FT_EEG.fsample = srate;
    
    % convert back to an EEGLAB file
    EEG = ft2eeglab(DETRENDED_FT_EEG);
        
    % now that detrending is complete, interpolate faulty electrodes
    if ~isempty(rejected_electrodes)
        EEG = pop_interp(EEG, orig_chanlocs, 'spherical');
        fid = fopen([outpath filesep 'bad_channels_' fname '.txt'], 'wt' );
        for c = 1:numel(rejected_electrodes)
            fprintf( fid, '%s\n', rejected_electrodes{c});
        end
        fclose(fid);
    end
    
    % can save as an FT_EEG file
    % save([outpath filesep fname '.mat'],'DETRENDED_FT_EEG');
    
    % better save dataset as EEGLAB file
    pop_saveset(EEG, 'filename',[fname '.set'],'filepath',outpath);
end


