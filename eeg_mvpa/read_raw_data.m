function [ FT_EEG, filename, chanlocs ]= read_raw_data(filepath,filename,outpath,channelset,erp_baseline,resample_eeg,do_csd,clean_data,clean_window,shuffle_trials,detrend_eeg)
% function [ FT_EEG, filename, chanlocs ]= read_raw_data(filepath,filename,outpath,channelset,erp_baseline,resample_eeg,do_csd,clean_data,clean_window,shuffle_trials,detrend_eeg)
% read eeglab data and transform into FT_EEG
% perform SCD and a bunch of other stuff such as trial rejection and resampling if necessary
% Johannes Fahrenfort, VU 2016, 2017
if nargin<11
    detrend_eeg = false;
end
if nargin<10
    shuffle_trials = false;
end
if nargin<9
    clean_window = [];
end
if nargin<8
    clean_data = false;
end
if nargin<7
    do_csd = false;
end
if nargin<6
    resample_eeg = false;
end
if nargin<5
    erp_baseline = 'no';
end
if nargin<4
    channelset = 'all';
end
if nargin<3
    outpath = filepath;
end
% remove .set if present
[~,filename,~] = fileparts(filename);


% first try to load eeglab data in .set format
if exist(fullfile(filepath,[filename '.set']),'file')
    EEG = pop_loadset('filename',[filename '.set'],'filepath',filepath);
    FT_EEG = eeglab2ft(EEG); % convert to internally used FT_EEG format
    chanlocs = EEG.chanlocs;
    clear EEG;
    % next, attempt to load fieldtrip data in .mat format in case this fails
elseif exist(fullfile(filepath,[filename '.mat']),'file')
    FT_EEG = load(fullfile(filepath,[filename '.mat']));
    FT_FIELDS = fieldnames(FT_EEG); % convert to internally used FT_EEG format
    if numel(FT_FIELDS) == 1
        FT_EEG = FT_EEG.(FT_FIELDS{1});
    end
else
    wraptext('The data needs to be stored in eeglab native format (.set extension) or in fieldtrip native format (.mat extension). The toolbox attempts to read your datafile using both extensions. If you are reading this message, both attempts have failed.');
    error(['Cannot load data, filename ''' filename ''' cannot be found at ''' filepath '''']);
end

% add sample info if missing
if ~isfield(FT_EEG,'fsample')
    FT_EEG.fsample = round((numel(FT_EEG.time)-1)/(FT_EEG.time(end)-FT_EEG.time(1)));
end

% double check whether all necessary fields are there
FT_FIELDS = fieldnames(FT_EEG);
reqfields = {'fsample', 'label', 'trial', 'time', 'trialinfo'};
absentfields = reqfields(~ismember(reqfields,FT_FIELDS));
if numel(absentfields)==1
    error(['The following required field is missing from your data: ' sprintf('%s',absentfields{1})]);
elseif numel(absentfields)>1
    error(['The following required fields are missing from your data: ' sprintf('%s, ',absentfields{1:end-1}) absentfields{end}]);
end

% detrend eeg
if detrend_eeg
    FT_EEG = detrend_FT_EEG(FT_EEG);
end

% baseline correction
cfg = [];
cfg.baseline = erp_baseline;
cfg.channel = 'all';
cfg.parameter = 'trial';
% turn off annoying FT warnings
if exist('ft_warning','file') == 2
    ft_warning off;
end
FT_EEG = ft_timelockbaseline(cfg,FT_EEG);
if clean_data % muscle artifact detection
    FT_EEG = muscle_removal(FT_EEG,[outpath filesep filename],clean_window);
end

% compute Current Source Density?
if do_csd
    cfg = [];
    cfg.channel = 'EEG';
    % turn off annoying FT warnings
    if exist('ft_warning','file') == 2
        ft_warning off;
    end
    FT_EEG = ft_selectdata(cfg,FT_EEG);
    % fix elec, which ft_selectdata should have done
    if isfield(FT_EEG,'elec')
        indx2rmv = ~ismember(FT_EEG.elec.label,FT_EEG.label);
        FT_EEG.elec.chanpos(indx2rmv,:) = [];
        FT_EEG.elec.elecpos(indx2rmv,:) = [];
        FT_EEG.elec.label(indx2rmv) = [];
    end
    cfg = [];
    cfg.method = 'spline'; % finite
    cfg.trials = 'all';
    FT_EEG = ft_scalpcurrentdensity(cfg,FT_EEG);
end

% resample EEG (this is only done in classify_RAW_eeglab_data, not in classify_TFR_from_eeglab_data)
if resample_eeg
    cfg = [];
    cfg.resamplefs = resample_eeg;
    % turn off annoying FT warnings
    if exist('ft_warning','file') == 2
        ft_warning off;
    end
    FT_EEG = ft_resampledata(cfg,FT_EEG);
end

% custom function to select channels
[~, channels] = select_channels(FT_EEG.label,channelset);
disp(['We had ' num2str(numel(FT_EEG.label)) ' channels/electrodes.']);
cfg = [];
cfg.channel = channels;
% turn off annoying FT warnings
if exist('ft_warning','file') == 2
    ft_warning off;
end
FT_EEG = ft_selectdata(cfg,FT_EEG);
disp(['Now we have ' num2str(numel(FT_EEG.label)) ' channels/electrodes.']);

% retain eeglab chanlocs for topoplots later on
if exist('chanlocs','var') % keep what came from eeglab
    [~, ~, chanindex] = intersect(FT_EEG.label,{chanlocs(:).labels},'stable'); % takes FT_EEG.label as point of departure!
    chanlocs = chanlocs(chanindex);
else % if not coming from eeglab, recreate eeglab chanlocs structure
    try
        chanfields = FT_FIELDS(ismember(FT_FIELDS,{'elec','grad'}));
        % if FT_EEG already contains chanloc data, use those data
        if ~isempty(chanfields)
            % remove field if not informative
            index2remove = [];
            for c = 1:numel(chanfields)
                if sum(ismember(FT_EEG.label,FT_EEG.(chanfields{c}).label)) == 0
                    FT_EEG = rmfield(FT_EEG,chanfields{c});
                    index2remove = [index2remove c];
                end
            end
            chanfields(index2remove) = [];
            % append all fields that are informative
            allchans = []; allpos = [];
            for c = 1:numel(chanfields)
                allchans = [allchans; FT_EEG.(chanfields{c}).label];
                allpos = [allpos; FT_EEG.(chanfields{c}).chanpos];
            end
            % restrict only to fields that have data and for which pos info exists
            for c = 1:numel(FT_EEG.label)
                findlabel = ismember(allchans,FT_EEG.label{c});
                if ~isempty(findlabel)
                    chanlocs(c).labels = FT_EEG.label{c};
                    chanlocs(c).X = allpos(findlabel,1);
                    chanlocs(c).Y = allpos(findlabel,2);
                    chanlocs(c).Z = allpos(findlabel,3);
                end
            end
        else
            % if no chanlocdata exist, read them in
            %chanlocdata = readlocs('plotting_1005.sfp','importmode','native');
            chanlocdata = readlocs(findcapfile,'importmode','native');
            [~, ~, chanindex] = intersect(FT_EEG.label,{chanlocdata(:).labels},'stable'); % takes FT_EEG.label as point of departure!
            chanlocs = chanlocdata(chanindex);
        end
    end
    if ~exist('chanlocs','var')
        disp('WARNING: could not find electrode or channel positions.');
        chanlocs = [];
    end
end

% destroy temporal order of trials
if shuffle_trials
    shuffle_order = randperm(numel(FT_EEG.trialinfo));
    FT_EEG.trialinfo = FT_EEG.trialinfo(shuffle_order);
    FT_EEG.trial = FT_EEG.trial(shuffle_order,:,:);
    FT_EEG.origindex = shuffle_order;
else
    FT_EEG.origindex = 1:numel(FT_EEG.trialinfo);
end