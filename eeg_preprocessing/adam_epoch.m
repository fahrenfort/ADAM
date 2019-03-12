function adam_epoch(cfg)
% Reads in CONTINUOUS (!) EEGLAB data, epochs and writes out as EPOCHED (!) EEGLAB file. 
% No pre-processing is applied other than mean removal from every trial after epoching.
% The output of this function can serve as input for ADAM_MVPA_FIRSTLEVEL.
%
%       cfg.datadir           = string specifiying the directory where the input files are
%                               located;
%       cfg.outputdir         = string specifying the directory where the results should be
%                               saved. Choose an informative name.
%       cfg.filenames         = N by 1 cell array containing N strings, in which each string
%                               contains the base of a filename (e.g. cfg.filenames =
%                               {'subject1', 'subject2'};). Do not add a file extension to the
%                               file name, ADAM will automatically look for EEGLAB .set files.
%       cfg.start_epoch       = in seconds (relative to target events)
%       cfg.end_epoch         = in seconds (relative to target events)
%       cfg.event_codes       = numeric array or comma separated list containing the relevant event
%                               values for the experiment. 
%
% Part of the ADAM toolbox, J.J.Fahrenfort, UvA/VU 2019

v2struct(cfg);

% do some checking
if ~exist('datadir','var')
    error('You need to specify in cfg.datadir where the data are located');
end
if ~exist('outputdir','var')
    error('You need to specify in cfg.outputdir where the results should be stored');
end
if ~exist('filenames','var')
    error('You need to specify a cell array in cfg.filenames containing the filenames containing the data of each subject');
end
if ~exist('event_codes','var')
    error('You need to specify in cfg.event_codes which event values you want to epoch on');
end
if ~exist('start_epoch','var')
    error('You need to specify in cfg.start_epoch defining the onset time of your epoch (in seconds)');
end
if ~exist('end_epoch','var')
    error('You need to specify in cfg.end_epoch defining the offset time of your epoch (in seconds)');
end

% some further input checking
if isempty(outputdir)
    outputdir = filepath;
else
    if ~exist(outputdir,'dir')
        mkdir(outputdir);
    end
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
if ischar(event_codes)
    conditions = regexp(event_codes, ',', 'split');
end
if isnumeric(event_codes)
    conditions = num2cell(event_codes);
end
% go
for filename = filenames
    [~,fname,~] = fileparts(filename{1});
    % load and filter
    EEG = pop_loadset('filename',[fname '.set'],'filepath',datadir);
    % epoch data
    if ~isempty(conditions)
        EEG = pop_epoch(EEG, conditions, [start_epoch end_epoch], 'newname', ['epoched_ ' fname], 'epochinfo', 'yes'); 
    end
    % remove mean of every trial
    EEG = pop_rmbase( EEG,[]);
    % save data
    pop_saveset(EEG, 'filename',[fname '.set'],'filepath',outputdir);
end


