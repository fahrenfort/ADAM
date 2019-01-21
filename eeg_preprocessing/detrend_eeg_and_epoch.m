function detrend_eeg_and_epoch(filepath,filenames,outpath, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, varargin)
% function detrend_eeg_and_epoch(filepath,filenames,outpath, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, varargin)
% reads in EEGLAB data, detrends, epochs and writes out. 
% filenames: names of files either in a cell array or as comma separated string. Wildcards * and ?
% can be used, e.g. filenames = '*' will take all the .set files in the input filepath as sources.
%
%       start_epoch             in seconds (relative to target events)
%       end_epoch               in seconds (relative to target events)
%       polynomial order        order of the polynomial that is used for detrending (Default: 30).
%       pad_length              wide epoch padding window used to fit the polynomial on, e.g.
%                               pad_length of 100 creates 50 sec pads around both sides of a trial
%                               (Default: 200). This padding window is only used during detrending
%                               and is discarded after the trial itsel is epoched out (when doing
%                               the final epoching).
%       start_mask              in seconds, specifies when your 'cognitive' or other events of
%                               interest start
%       end_mask                in seconds, specifies when your 'cognitive' or other events of
%                               interest end 
%       mask_only_current       false or true (default: false). Specifies whether the algorithm for
%                               any given trial masks out all other trials from the wide window pad,
%                               or only the currently relevant trial.
% 
% varargin is list of conditions to epoch, either as strings or numbers
% separated by comma's, or as a cell array
% Use  cond_string([1001:1024 1101:1124]) when called from create_qsub_files in
% order to epoch all conditions running from 1001 to 1024 and from 1101 to
% 1124. 
%
% example: detrend_eeg_and_epoch('c:\alldata\','subject01', 'c:\filtereddata', -.5, 1, 30, 100, 0, .6, cond_string([1001:1024 1101:1124]));
%
% J.J.Fahrenfort, VU 2018, 2019

% some input checking
filepath,filenames,outpath, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, 

if nargin < 10
    mask_only_current = false;
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
    pad_length = 200;
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
    if ~iscell(varargin{1})
        conditions = regexp(varargin{1}, ',', 'split');
    end
else
    conditions = varargin;
    if ~isempty(conditions)
        if ~ischar(conditions{1})
            conditions = regexp(vec2str([conditions{:}]), ',', 'split');
        end
    end
end

% go
for filename = filenames
    [~,fname,~] = fileparts(filename{1});
    
    % load EEGLAB data and immediately convert to FT_EEG format for ease of processing
    FT_EEG = eeglab2ft(fname,filepath,true); % third argument indicated data is continuous
    
    % obtain names to work with
    eeg_data = FT_EEG.trial{1};
    eeg_time = FT_EEG.time{1};
    srate = FT_EEG.fsample;
    trialinfo = FT_EEG.trialinfo;
    label = FT_EEG.label;
    
    % select relevant events
    trialinfo = trialinfo(ismember(trialinfo(:,1),conditions),:);
    
    % mirror-pad edges of the unepoched data, so that extracting wide padded epochs will not be problematic
    eeg_data = padarray(eeg_data,[0 pad_length*srate],'both','symmetric');
    eeg_time = padarray(eeg_time,[0 pad_length*srate],NaN,'both');
        
    % create a mask for all trials
    eeg_mask = ones(size(eeg_data));
    if ~mask_only_current
        for cTrials = 1:size(trialinfo,1)
            mask_startind = nearest(eeg_time,trialinfo(cTrials,2)+1000*start_mask);
            mask_stopind = nearest(eeg_time,trialinfo(cTrials,2)+1000*stop_mask);
            eeg_mask(mask_startind:mask_stopind) = 0;
        end
    end
    
    for cTrials = 1:size(trialinfo,1)
        
        % identify zero point and borders of trial
        start_ind = nearest(eeg_time,trialinfo(cTrials,2)+start_epoch*1000);
        zero_ind = nearest(eeg_time,trialinfo(cTrials,2));
        stop_ind = nearest(eeg_time,trialinfo(cTrials,2)+end_epoch*1000);
        
        % extract wide padded trial
        pad_time = eeg_time((start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
        pad_data = eeg_data((start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
        if mask_only_current % mask only current trial
            pad_mask = ones(size(pad_data));
            mask_zero_ind = nearest(pad_time,trialinfo(cTrials,2)+1000*start_mask);
            mask_stop_ind = nearest(pad_time,trialinfo(cTrials,2)+1000*stop_mask);
            pad_mask(mask_zero_ind:mask_stop_ind) = 0;
        else % mask all trials
            pad_mask = eeg_mask((start_ind-pad_length/2*srate):(stop_ind+pad_length/2*srate));
        end
        
        % create a mask matrix
        wt = repmat(pad_mask,[numel(label) 1]);
        
        % new: make sure data is physically removed before calling nt_detrending function
        pad_data_orig = pad_data;
        
        % no way in hell that the detrending function can make use of brain data:
        pad_data_nocogdata = pad_data .* wt;
        
        % call nt_detrend without brain data
        [tmp,w1,~,regressline] = nt_detrend(pad_data_nocogdata',1,wt'); % start with 1st order
        
        % take out the regression line manually
        pad_data = pad_data_orig - regressline';
        
        % call nt_detrend with higher order polynomial
        [x2, w2,~,regressline] = nt_detrend(tmp,polynom,w1); % then nth order with mask of previous step
        
        % take out the regression line manually
        clean_data = pad_data-regressline'; % manually take out regression slope from actual data
        
        % go back to narrow epoch
        trial(cTrial,:,:) = clean_data(:,(pad_length/2*srate+1):end-pad_length/2*srate);
        new_time = pad_time - eeg_time(zero_ind);
        new_time = new_time((pad_length/2*srate+1):end-pad_length/2*srate);
    end
    
        % create a new epoched FT_EEG dataset
%     dimord: 'rpt_chan_time'
%     label: {27x1 cell}
%     origindex: [1x574 double]
%     time: [1x750 double]
%     trial: [574x27x750 double]
%     trialinfo: [574x1 double]
    NEW_FT_EEG.dimord = 'rpt_chan_time';
    NEW_FT_EEG.label = label;
    NEW_FT_EEG.time = new_time;
    NEW_FT_EEG.trial = trial;
    NEW_FT_EEG.trialinfo = trialinfo(:,1);
    
    % create and save a figure with:
    % - butterfly ERP plot across epoch window prior to detrending
    % - butterfly ERP plot across epoch window after detrending
    % - single trial halfway the experiment containing butterfly plot prior to detrending
    % - single trial halfway the experiment containing butterfly plot after detrending
    
    % save dataset
    save([outpath filesep fname '.mat'],'NEW_FT_EEG');
end


