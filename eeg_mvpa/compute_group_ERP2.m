function group_FT_EEG = compute_group_ERP2(folder_name,cfg)
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
        group_FT_EEG{cdirz} = subcompute_group_ERP([folder_name filesep dirz{cdirz}],cfg);
    end
    cfg.folder = folder_name;
else
    if ~exist('folder_name','dir') && ~iscell(folder_name) 
        error([folder_name ' should refer to a full and existing folder path. Alternatively leave folder_name empty to pop up a selection dialog.']);
    end
    group_FT_EEG  = subcompute_group_ERP(folder_name,cfg);
end

% fill cfg.plot_order in case not given by user and from single folder
if numel(cfg.plot_order) == 1
    cfg.plot_order = {stats(:).condname};
end

function group_FT_EEG = subcompute_group_ERP(folder_name,cfg)
channelpool = 'ALL';
% unpack graphsettings
v2struct(cfg);

% get filenames
plotFreq = ''; % this is empty for now, but might be used to look at the ERPs obtained from a TF analysis
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

% do the loop, restrict time and frequency if applicable
for cSubj = 1:nSubj
    fprintf(1,'loading subject %d of %d\n', cSubj, nSubj);
    matObj = matfile([folder_name filesep channelpool plotFreq filesep subjectfiles{cSubj}]);
    if ~isempty(whos(matObj,'settings'))
        settings = matObj.settings;
    else % BW compatibility, will become obsolete over time
        error('Cannot find settings struct which the contains the ERP info');
    end
    % OLD for backwards compatible
    if isfield(settings,'ft_erpstruct')
        FT_ERP = settings.ft_erpstruct;
    elseif isfield(settings,'FT_ERP')
        FT_ERP = settings.FT_ERP;
    else
        error('no ERP was computed for these data');
    end
    
    if iscell(FT_ERP)
        FT_ERP = FT_ERP{2};
        disp('Extracting ERPs from testing data');
    end
    group_FT_EEG{cSubj} = FT_ERP;
    
end


    
