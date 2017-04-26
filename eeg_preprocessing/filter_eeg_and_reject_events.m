function filter_eeg_and_reject_events(filepath,filenames,outpath, locutoff, hicutoff, start_epoch, end_epoch, varargin)
% function filter_eeg_and_epoch(filepath,filename,outpath, locutoff, hicutoff, start_epoch, end_epoch, varargin)
% reads in EEGLAB data, filters, and writes out.
% The portions of the data that contain events are REMOVED, so that the
% resulting data can be used as a basis for ICA blink detection without
% running the risk of having cognitive events (i.e. ERPs) feeding into
% blink components.
% filenames: names of files either in a cell array or as comma separated
% string. Wildcards * and ? can be used, e.g. filenames = '*' will
% take all the .set files in the input filepath as sources. 
%
%     locutoff  - lower edge of the frequency pass band (Hz)
%                 {[]/0 -> lowpass} 
%     hicutoff  - higher edge of the frequency pass band (Hz)
%                 {[]/0 -> highpass}
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
% J.J.Fahrenfort, VU 2015

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
if ~isempty(strfind(filenames,'*')) || ~isempty(strfind(filenames,'?'))
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
    % load and filter
    EEG = pop_loadset('filename',[fname '.set'],'filepath',filepath);
    if any([locutoff hicutoff]>0)
        EEG = pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder, revfilt);
    end
    % remove mean across the entire block (somewhat superfluous but can't hurt)
    EEG = pop_rmbase( EEG,[]);
    % epoch data and parts of the continuous data that fall within epochs
    if ~isempty(conditions)
        latencies2reject = [];
        for cEvent = 1:numel(EEG.event)
            if ~ischar(EEG.event(cEvent).type)
                type = num2str(EEG.event(cEvent).type);
            else
                type = EEG.event(cEvent).type;
            end
            if any(strcmp(type,conditions))
                start = EEG.event(cEvent).latency + start_epoch*1000;
                stop = EEG.event(cEvent).latency + end_epoch*1000;
                if start > 0 && stop < EEG.pnts
                    latencies2reject = [latencies2reject; start stop];
                end
            end
        end
        EEG = eeg_eegrej( EEG, latencies2reject);
    end
    % save data
    pop_saveset(EEG, 'filename',[fname '.set'],'filepath',outpath);
end


