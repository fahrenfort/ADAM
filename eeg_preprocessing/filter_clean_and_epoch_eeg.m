function filter_clean_and_epoch_eeg(filepath,filenames,outpath, locutoff, hicutoff, start_epoch, end_epoch, varargin)
% function filter_eeg_and_epoch(filepath,filename,outpath, locutoff, hicutoff, start_epoch, end_epoch, varargin)
% reads in EEGLAB data, filters, removes and interpolates bad channels, epochs, and writes out.
% filenames: names of files either in a cell array or as comma separated
% string. Wildcards * and ? can be used, e.g. filenames = '*' will
% take all the .set files in the input filepath as sources.
%
%     locutoff  - lower edge of the frequency pass band (Hz)
%                 {[]/0 -> only lowpass}
%     hicutoff  - higher edge of the frequency pass band (Hz)
%                 {[]/0 -> only highpass}
%
% if both locutoff and hicutoff are zero the function does not filter
% varargin is list of conditions to epoch, either as strings or numbers
% separated by comma's, or as a cell array
% Use  vec2str([1001:1024 1101:1124]) when called from create_qsub_files in
% order to epoch all conditions running from 1001 to 1024 and from 1101 to
% 1124. If varargin is left empty, no epoching will be applied.
%
% example: filter_eeg_and_epoch('c:\alldata\','subject01', 'c:\filtereddata', 0.1, 0, -.25, 0, num2str([1001:1024 1101:1124]));
%
% J.J.Fahrenfort, VU 2015, UvA 2017

% some input checking
revfilt = 0;
filtorder = [];
if nargin < 5
    hicutoff = 0;
end
if nargin < 4
    locutoff = 0;
end
if nargin < 3
    error('not enough arguments');
end
if isempty(outpath)
    outpath = filepath;
else
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
end
if ischar(hicutoff);
    hicutoff = string2double(hicutoff);
end
if ischar(locutoff)
    locutoff = string2double(locutoff);
end
if ischar(start_epoch);
    start_epoch = string2double(start_epoch);
end
if ischar(end_epoch)
    end_epoch = string2double(end_epoch);
end
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
    % load 
    EEG = pop_loadset('filename',[fname '.set'],'filepath',filepath);
    % hp filter
    if locutoff>0
        EEG = pop_eegfiltnew(EEG, locutoff, [], filtorder, revfilt);
    end
    % find EEG channels
    channelnames = {EEG.chanlocs(:).labels};
    try
        eog_channels = select_channels(channelnames,'EOG');
        eog_present = true;
    catch
        disp('No EOG channels in data');
        eog_present = false;
    end
    try 
        eeg_channels = select_channels(channelnames,'EEG');
    catch
        error('stopping now, there are no EEG channels in this set??');
    end
    % double check whether channel location information is present
    nopos_channels = [];
    for cEl=1:length(EEG.chanlocs)
        if (any(isempty(EEG.chanlocs(1,cEl).X)&isempty(EEG.chanlocs(1,cEl).Y)&isempty(EEG.chanlocs(1,cEl).Z)&isempty(EEG.chanlocs(1,cEl).theta)&isempty(EEG.chanlocs(1,cEl).radius)))
            nopos_channels = [nopos_channels cEl];
        end
    end
    % separate EEG and EOG channels
    if eog_present
        EOG = pop_select(EEG, 'channel', eog_channels);
    end
    EEG = pop_select(EEG, 'channel', eeg_channels);
    % look up electrode info
    if any(ismember(eeg_channels,nopos_channels))
        disp(['WARNING: Channels ' num2str(nopos_channels) ' have incomplete location information. Now attempting to read in location information.']);
        EEG = pop_chanedit(EEG, 'lookup', trycapfile);
    end
    % identify and remove bad channels from the EEG, interpolate and write bad channels to text file
    chanlocs = EEG.chanlocs;
    EEG = clean_artifacts(EEG,'Highpass','off','WindowCriterion','off','burst_crit','off');
    realchanlocs = EEG.chanlocs;
    rejected_electrodes = setdiff({chanlocs.labels},{realchanlocs.labels});
    if ~isempty(rejected_electrodes)
        EEG = pop_interp(EEG, chanlocs, 'spherical');
        fid = fopen([outpath filesep 'bad_channels_' fname '.txt'], 'wt' );
        for c = 1:numel(rejected_electrodes)
            fprintf( fid, '%s\n', rejected_electrodes{c});
        end
        fclose(fid);
    end
    % put EOG channels back in
    if eog_present
        EEG.nbchan = EEG.nbchan + EOG.nbchan;
        EEG.data(end+1:end+EOG.nbchan,:,:) = EOG.data;
        EEG.chanlocs(end+1:end+EOG.nbchan) = EOG.chanlocs;
    end
    EEG = eeg_checkset(EEG);
    % lowpass
    if hicutoff>0
        EEG = pop_eegfiltnew(EEG, [], hicutoff, filtorder, revfilt);
    end
    % epoch data
    if ~isempty(conditions)
        EEG = pop_epoch(EEG, conditions, [start_epoch end_epoch], 'newname', ['epoched_ ' fname], 'epochinfo', 'yes');
    end
    % remove mean across the epochs (somewhat superfluous but can't hurt)
    EEG = pop_rmbase( EEG,[]);
    % save data
    pop_saveset(EEG, 'filename',[fname '.set'],'filepath',outpath);
end


