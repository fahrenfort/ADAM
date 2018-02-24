function [FT_EEG, boollist] = muscle_removal(EEG,outfile,clean_window,cutoff)
% [FT_EEG, boollist] = muscle_removal(EEG,outfile,clean_window,cutoff)
% does muscle artefact detection using fieldtrip
% takes FT_EEG or EEGLAB as input format
% outfile is the file (including path) where the image containing the
% output should be saved
% clean_window is in seconds
% cutoff is the z-value to use as cutoff, uses -3 by default (see ft_artifact_zvalue_jjf)
% outputs muscle_boolean of trials to be removed (1 = remove) as well as a
% cleaned FT_EEG struct
% boollist contains 1 for the trials that have a muscle artefact
%
% J.J.Fahrenfort, VU 2016

% NEED TO FIGURE OUT HOW cfg.padtype = 'mirror' WORKS, SEE http://www.fieldtriptoolbox.org/reference/ft_preprocessing

if nargin<4
    cutoff = -3;
end
if nargin<3
    clean_window = [];
end
if nargin<2
    outfile = [];
end
if ~isempty(outfile)
    [filepath,filename,~] = fileparts(outfile);
    outfile = fullfile(filepath,['muscle_artefacts_' filename]);
end

if isfield(EEG,'setname')
    FT_EEG = eeglab2ft(EEG);
else
    FT_EEG = EEG;
end
clear EEG;

try
    % create sample info (no actual baseline correction is performed)
    cfg = [];
    cfg.channel = 'all';
    cfg.parameter = 'trial';
    FT_EEG = ft_timelockbaseline(cfg,FT_EEG);
    
    % compute fsample
    if ~exist('fsample','var')
        fsample = (numel(FT_EEG.time)-1)/(FT_EEG.time(end) - FT_EEG.time(1));
    end
    
    % pick the window that should be cleaned, and add some padding around it
    if ~isempty(clean_window)
        cfg = [];
        padding = min([clean_window(1) - min(FT_EEG.time), max(FT_EEG.time) - clean_window(2)]); % take the maximum padding allowed
        if padding < 0
            error('cannot clean data without some padding, please include a baseline window in your data prior to you window of interest');
        end
        cfg.latency = [clean_window(1) - padding, clean_window(2) + padding];
        FT_TEMP = ft_selectdata(cfg,FT_EEG);
    else
        FT_TEMP = FT_EEG; % if no window is given, pad the baseline (pre-0) on both sides
        padding = -min(FT_EEG.time);
    end
    
    % cutoff and padding, only use certain electrodes
    cfg             = [];
    cfg.continuous  = 'no';
    % select only EEG electrodes
    % [~, electrodes] = select_channels(FT_EEG.label,'EEG');
    cfg.artfctdef.zvalue.channel    = FT_EEG.label; % all channels
    cfg.artfctdef.zvalue.cutoff     = cutoff; % was 15
    cfg.artfctdef.zvalue.trlpadding = -padding; % exclude parts of the trial on both ends
    cfg.artfctdef.zvalue.fltpadding = 0; % this setting does not always work with every bdf, then set to 0
    cfg.artfctdef.zvalue.artpadding = 0; %
    % algorithmic parameters
    % cfg.artfctdef.zvalue.method     = 'trialdemean'; % stop doing this
    % cfg.artfctdef.zvalue.polyremoval = 'yes';
    cfg.artfctdef.zvalue.bpfilter   = 'yes';
    % make sure that frequency range is reasonable given nyquist
    freqrange = [110 140];
    if freqrange(1)/fsample > .5 % remember nyquist
        disp('WARNING: cannot do muscle rejection using this sampling rate');
        boollist = zeros(size(FT_EEG.trialinfo));
        return
    end
    if freqrange(2)/fsample > .5 % remember nyquist
        disp('WARNING: lowering max frequency that will be searched for muscle detection');
        freqrange(2) = floor(fsample/2)-1;
    end
    cfg.artfctdef.zvalue.bpfreq     = freqrange;
    cfg.artfctdef.zvalue.bpfiltord  = 6; % was 9, if it doesn't work try 6
    cfg.artfctdef.zvalue.bpfilttype = 'but';
    cfg.artfctdef.zvalue.hilbert    = 'yes';
    cfg.artfctdef.zvalue.boxcar     = 0.2;
    % make the process interactive?
    cfg.artfctdef.zvalue.interactive= 'no';
    if ~isempty(outfile)
        cfg.saveoutput = outfile;
    end
    % go
    [cfg, artifact]  = ft_artifact_zvalue_jjf(cfg,FT_TEMP);
    clear FT_TEMP;
    % cfg = ft_databrowser(cfg,FT_EEG); % turn this on to inspect artefacts manually
    % generate a list of artifact indices
    smpl = cfg.artfctdef.zvalue.trl(:,1:2);
    boollist = zeros(size(smpl,1),1);
    for c=1:size(artifact,1)
        boollist1 = artifact(c,1) > smpl(:,1) & artifact(c,1) < smpl(:,2);
        boollist2 = artifact(c,2) > smpl(:,1) & artifact(c,2) < smpl(:,2);
        boollist = boollist | boollist1 | boollist2;
    end
    % remove artefacts
    FT_EEG.trialinfo = FT_EEG.trialinfo(~boollist,:);
    if isfield(FT_EEG,'sampleinfo')
        FT_EEG.sampleinfo = FT_EEG.sampleinfo(~boollist,:);
    end
    FT_EEG.trial = FT_EEG.trial(~boollist,:,:);
    % keep track
    trialindex = find(boollist);
    if ~isempty(outfile)
        fid = fopen([outfile '.txt'], 'wt' );
        for c = 1:numel(trialindex)
            fprintf( fid, '%d\n', trialindex(c));
        end
        fclose(fid);
    end
    boollist = boollist';
catch ME
    disp('WARNING: could not remove muscle artefacts due to an error:');
    disp(ME);
    boollist = [];
end