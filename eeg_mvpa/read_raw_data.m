function [ FT_EEG, filename, chanlocs ]= read_raw_data(filepath,filename,outpath,msettings)
% function [ FT_EEG, filename, chanlocs ]= read_raw_data(filepath,filename,outpath,msettings)
% read eeglab data or FT data and transform into FT_EEG
% perform CSD and a bunch of other stuff such as trial rejection and resampling if necessary
% internal function of the ADAM toolbox.
% Johannes Fahrenfort, VU 2016, 2017, 2018

% set some defaults
shuffle_trials = false;
clean_window = [];
clean_data = false;
do_csd = false;
resample_eeg = false;
resample_method = 'resample';
erp_baseline = 'no';
channelpool = 'all';
if nargin>3
    v2struct(msettings);
end
if nargin<3 || isempty(outpath)
    outpath = filepath;
end

% remove .set / .mat if present
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
    % IMIPLEMENT CONVERSION TO trial x electrode x time
else
    wraptext('The data needs to be stored in eeglab native format (.set extension) or in fieldtrip native format (.mat extension). The toolbox attempts to read your datafile using both extensions. If you are reading this message, both attempts have failed.');
    error(['Cannot load data, filename ''' filename ''' cannot be found at ''' filepath '''']);
end

% convert to FT timelock
cfg = [];
cfg.baseline = [];
cfg.channel = 'all';
cfg.parameter = 'trial';
% turn off annoying FT warnings
if exist('ft_warning','file') == 2
    ft_warning off;
end
FT_EEG = ft_timelockbaseline(cfg,FT_EEG);

% baseline correction, converts data to FT timelock format
cfg = [];
cfg.baseline = erp_baseline;
cfg.channel = 'all';
cfg.parameter = 'trial';
% turn off annoying FT warnings
if exist('ft_warning','file') == 2
    ft_warning off;
end
FT_EEG = ft_timelockbaseline(cfg,FT_EEG);

% add dimord info if missing
if ~isfield(FT_EEG,'dimord')
    if size(FT_EEG.trial,1) == numel(FT_EEG.trialinfo) && size(FT_EEG.trial,2) == numel(FT_EEG.label) && size(FT_EEG.trial,3) == numel(FT_EEG.time)
        FT_EEG.dimord = 'rpt_chan_time';
    else
        error('The dimord field seems missing and/or dimensions in the dataset seem to be off. When you are using fieldtrip format, make sure the input files have dimord rpt_chan_time (trial x channel x time).');
    end
end

% keep track of index numbers
origindex = 1:numel(FT_EEG.trialinfo);
% muscle artifact detection
if clean_data 
    [FT_EEG, boollist] = muscle_removal(FT_EEG,[outpath filesep filename],clean_window);
    origindex = origindex(~boollist);
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
    fsample = (numel(FT_EEG.time)-1)/(FT_EEG.time(end) - FT_EEG.time(1));
    if mod(fsample,resample_eeg) ~=0 && ~strcmpi(resample_method,'resample')
        disp('Cannot use downsample or average_timebin method because new sampling rate is not a proper divisor of original sampling rate. Reverting to downsample method.')
        resample_method = 'resample';
    end
    if strcmpi(resample_method,'average_timebin')
        disp(['Original sampling rate is ' num2str(fsample) ' Hz.']);
        disp('Resampling by averaging time samples in bin (no anti-alias filter).');
        nsamples = fsample/resample_eeg;
        centersample = ceil((nsamples+1)/2); % conservative estimate of center (always averaging across the past more than across the future)
        minindex = centersample - 1;
        plusindex = nsamples - centersample;
        newtime = FT_EEG.time(centersample:nsamples:end);
        cnew = 0;
        newtrial= zeros(size(FT_EEG.trial,1),size(FT_EEG.trial,2),numel(newtime));
        for csample = centersample:nsamples:size(FT_EEG.trial,3)
            cnew = cnew + 1;
            newtrial(:,:,cnew) = mean(FT_EEG.trial(:,:,(csample-minindex):(csample+plusindex)),3);
        end
        FT_EEG.time = newtime;
        FT_EEG.trial = newtrial;
        fsample = (numel(FT_EEG.time)-1)/(FT_EEG.time(end) - FT_EEG.time(1));
        disp(['New sampling rate is ' num2str(fsample) ' Hz.']);
    else
        cfg = [];
        cfg.resamplefs = resample_eeg;
        cfg.resamplemethod = resample_method;
        cfg.detrend = 'no';
        % turn off annoying FT warnings
        if exist('ft_warning','file') == 2
            ft_warning off;
        end
        FT_EEG = ft_resampledata(cfg,FT_EEG);
    end
end

% custom function to select channels
[~, channels] = select_channels(FT_EEG.label,channelpool);
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
            disp('WARNING: using electrode positions that are native to the data set. Therefore, the direction of the nose in topoplots cannot be ascertained with certainty. If needed, you can adjust the cfg.nosedir property prior to plotting (see help adam_plot_BDM_weights).');
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
            if exist('1005chanlocdata.mat','file')
                load('1005chanlocdata.mat');
            else
                chanlocdata = readlocs(trycapfile,'importmode','native'); % from standard 10-20 system
            end
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
    FT_EEG.origindex = origindex(shuffle_order);
else
    FT_EEG.origindex = origindex;
end

% add sample info if missing
if ~isfield(FT_EEG,'fsample')
    FT_EEG.fsample = round((numel(FT_EEG.time)-1)/(FT_EEG.time(end)-FT_EEG.time(1)));
end
% double check whether all necessary fields are now present
FT_FIELDS = fieldnames(FT_EEG);
reqfields = {'fsample', 'label', 'trial', 'time', 'trialinfo', 'dimord'};
absentfields = reqfields(~ismember(reqfields,FT_FIELDS));
if ~isempty(absentfields)
    error(['The following required fields are missing from your data: ' cell2csv(absentfields,true)]);
end

% remove superfluous fields
FT_FIELDS = fieldnames(FT_EEG);
reqfields = {'fsample', 'label', 'trial', 'time', 'trialinfo', 'dimord', 'origindex'};
rmfields = setdiff(FT_FIELDS,reqfields);
FT_EEG = orderfields(rmfield(FT_EEG,rmfields));