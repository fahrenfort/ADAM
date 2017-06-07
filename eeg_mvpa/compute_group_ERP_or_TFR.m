function group_FT_EEG = compute_group_ERP_or_TFR(folder_name,cfg)
% Get all the ERPs from settings and put them into a cell array for further
% analysis in fieldtrip
%
% By J.J.Fahrenfort, VU, 2017

% backwards compatibility
v2struct(cfg);
if exist('plotorder','var')
    plot_order = plotorder;
    cfg.plot_order = plot_order;
    cfg = rmfield(cfg,'plotorder');
end

% Main routine, is a folder name specified? If not, pop up selection dialog
if isempty(folder_name)
    if ~isfield(cfg,'startdir')
        cfg.startdir = '';
        disp('NOTE: it is easier to select a directory when you indicate a starting directory using cfg.startdir, otherwise you have to start selection from root every time...');
    end
    folder_name = uigetdir(cfg.startdir,'select directory to plot');
    if ~ischar(folder_name)
        return
    end
    if ~exist('plot_order','var') || isempty(plot_order)
        dirz = dir(folder_name);
        dirz = {dirz([dirz(:).isdir]).name};
        dirz = dirz(cellfun(@isempty,strfind(dirz,'.')));
        cfg.plot_order = dirz;
    else
        dirz = cfg.plot_order;
    end
    % loop through directories
    for cdirz = 1:numel(dirz)
        group_FT_EEG{cdirz} = subcompute_group([folder_name filesep dirz{cdirz}],cfg);
    end
    cfg.folder = folder_name;
else
    if ~exist('folder_name','dir') && ~iscell(folder_name) 
        error([folder_name ' should refer to a full and existing folder path. Alternatively leave folder_name empty to pop up a selection dialog.']);
    end
    group_FT_EEG  = subcompute_group(folder_name,cfg);
end

% fill cfg.plot_order in case not given by user and from single folder
if numel(cfg.plot_order) == 1
    cfg.plot_order = {stats(:).condname};
end

function group_FT_EEG = subcompute_group(folder_name,cfg)
channelpool = 'ALL';
% unpack graphsettings
v2struct(cfg);

% get filenames
freqfolder_contains_freq_time = ~isempty(dir([folder_name filesep channelpool filesep 'allfreqs']));
if freqfolder_contains_freq_time
    datafield = 'TFR';
    plotFreq = [filesep 'allfreqs'];
else
    datafield = 'ERP';
    plotFreq = '';
end
cfg.plotFreq = plotFreq;
subjectfiles = dir([folder_name filesep channelpool plotFreq filesep '*.mat']);
[~, condname] = fileparts(folder_name);
name = strsplit(folder_name,filesep);
name = name(1:end-1);
subjectfiles = { subjectfiles(:).name };
nSubj = numel(subjectfiles);
if nSubj == 0
    error(['cannot find data in specified folder ' folder_name filesep channelpool plotFreq filesep]);
end

% do the loop
for cSubj = 1:nSubj
    fprintf(1,'loading subject %d of %d\n', cSubj, nSubj);
    matObj = matfile([folder_name filesep channelpool plotFreq filesep subjectfiles{cSubj}]);
    
    % get goodies
    if strcmpi(datafield,'ERP')
        FT_EEG = matObj.FT_ERP;
    else
        FT_EEG = matObj.FT_TFR;
    end
    
    if iscell(FT_EEG)
        FT_EEG = FT_EEG{2};
        disp(['Extracting ' datafield 's from testing data']);
    end
    group_FT_EEG{cSubj} = FT_EEG;
    
end


    
