function [TFR, groupTFRhann, groupTFRmult] = compute_TFR(FT_EEG,sampling_rate,frequencies)
% function TFR = compute_TFR(FT_EEG,sampling_rate,frequencies)
% takes standard FT input and computes single trial TFR for MVPA and
% separate group TFRs (no single trial data) for Hanning and multitaper 
% FT_EEG is standard fieldtrip struct
% sampling_rate is the frequency to which the output should be downsampled
% (if left empty, output has the same sr as the input
% freqs is a specification of the frequencies to include
%

if nargin<3
    frequencies = 2:2:100;
end
if isempty(frequencies)
    frequencies = 2:2:100;
end

% downsample_factor specifies the factor by which to downsample
if nargin < 2
    sampling_rate = 0;
end
if ischar(sampling_rate)
    sampling_rate = string2double(sampling_rate);
end
if sampling_rate == 0
    sampling_rate = FT_EEG.fsample;
end
downsample_factor = FT_EEG.fsample/sampling_rate;
% check if downsample_factor is an integer
if ~(floor(downsample_factor)==downsample_factor)
    downsample_factor = round(downsample_factor);
    disp('WARNING: Specified resampling_rate resulted in non-integer downsample_factor: rounding off new sampling rate to prevent interpolation. Next time you better pick a sampling rate that fits.');
end
% times of interest, possibly better to keep actual time points than to create new time series
if iscell(FT_EEG.time)
    toi = downsample(FT_EEG.time{1},downsample_factor); 
else
    toi = downsample(FT_EEG.time,downsample_factor); 
end

% windows and cycles of interest (use appproximate time-window of .5 up to
% 30 Hz, above that, we use frequency-dependent time-window, by keeping the
% cycles constant at 15) 
% Few cycles: good temporal, poor frequency precision, Many cycles: bad temporal, good frequency precision
% low frequencies: few cycles, high frequencies: many cycles
for cFreq=1:numel(frequencies)
    for cCycle=1:100
        tim_winspace(cFreq,cCycle) = cCycle/frequencies(cFreq);
        % we use a time window of approx .5 seconds to keep temporal
        % resolution high in the lower frequencies
        % then from 30 Hz onwards we keep the nr of cycles constant to
        % increase temporal resolution
        act_cyclespace(cFreq) = nearest(tim_winspace(cFreq,:),.5); 
        act_ftimwin(cFreq) = tim_winspace(cFreq,act_cyclespace(cFreq));
    end
end

% get single condition data for group analysis
condlabels = unique(FT_EEG.trialinfo);
for c = 1:numel(condlabels)
    cfg              = [];
    cfg.trials       = FT_EEG.trialinfo == condlabels(c);
    FT_COND{c}       = ft_selectdata(cfg,FT_EEG);
end

% Hanning taper for low frequencies
cfg              = [];
cfg.output       = 'pow';
%cfg.channel      = 'EEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = frequencies(frequencies<=30);                        % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = act_ftimwin(frequencies<=30); % repmat(.5,numel(cfg.foi),1);   % length of time window = 0.5 sec
cfg.toi          = toi;                          % the times around which TFRs are computed, determined by downsample_factor
cfg.pad          = 'maxperlen';
if ~isempty(cfg.foi)
    % single trial
    cfg.keeptrials   = 'yes';
    TFRhann          = ft_freqanalysis(cfg, FT_EEG);
    TFRhann.cumtapcnt = TFRhann.cumtapcnt(1,:); % remove redundant information, number of tapers is the same for each trial
    % also get condition data on group level
    for c = 1:numel(condlabels)
        cfg.keeptrials   = 'no';
        groupTFRhann.(['cond_' num2str(condlabels(c))])  = ft_freqanalysis(cfg, FT_COND{c});
    end
else
    groupTFRhann     = [];
end
   
% multitaper for higher frequencies
cfg = [];
cfg.output     = 'pow';
% cfg.channel    = 'EEG';
cfg.method     = 'mtmconvol';
cfg.taper      = 'dpss';
cfg.foi        = frequencies(frequencies>30);
cfg.t_ftimwin  = 15./cfg.foi; % -> this decreases window size for higher freqs, keeping cycles constant at 15
cfg.tapsmofrq  = 0.1 * cfg.foi; % -> this increases frequency smoothing for higher freqs
%cfg.t_ftimwin  = repmat(.4,numel(cfg.foi),1); % -> keeps window size constant at .4 as in PNAS paper
%cfg.tapsmofrq  = 10; % -> keeps frequency smoothing constant at 20 Hz -> result: 7 tapers as in PNAS paper
cfg.toi        = toi;
cfg.pad        = 'maxperlen';
cfg.keeptrials = 'yes';
if ~isempty(cfg.foi)
    cfg.keeptrials   = 'yes';
    TFRmult        = ft_freqanalysis(cfg, FT_EEG);
    TFRmult.cumtapcnt = TFRmult.cumtapcnt(1,:); % remove redundant information, number of tapers is the same for each trial
    % also get condition data on group level
    for c = 1:numel(condlabels)
        cfg.keeptrials   = 'no';
        groupTFRmult.(['cond_' num2str(condlabels(c))])  = ft_freqanalysis(cfg, FT_COND{c});
    end
else
    groupTFRmult = [];
end

% free up memory
clear FT_EEG FT_COND;

% get relevant part of the spectrum
if exist('TFRhann','var')
    TFR             = TFRhann; 
else
    TFR             = TFRmult;
end
% if both exist, merge the results
if exist('TFRhann','var') && exist('TFRmult','var')
    TFR             = rmfield(TFR,'cfg'); % remove the cfg field
    TFR.freq        = cat(2,TFRhann.freq,TFRmult.freq);
    TFR.powspctrm   = cat(3,TFRhann.powspctrm,TFRmult.powspctrm);
    TFR.cumtapcnt   = cat(2,TFRhann.cumtapcnt,TFRmult.cumtapcnt);
    TFR.cfg.hann    = TFRhann.cfg;
    TFR.cfg.mult    = TFRmult.cfg;
end
% finally cut out annoying NaNs due to padding, a bit complex due to unknown format
dims = regexp(TFR.dimord, '_', 'split');
chandim = find(strcmp(dims,'chan'));
timedim = find(strcmp(dims,'time'));
trialdim = find(strcmp(dims,'rpt'));
freqdim = find(strcmp(dims,'freq'));
index    = cell(1, ndims(TFR.powspctrm));
index(:) = {':'};
index{chandim} = 1;
index{trialdim} = 1;
index{timedim} = find(~isnan(squeeze(sum(TFR.powspctrm(index{:}),freqdim))));
index{chandim} = ':';
index{trialdim} = ':';
TFR.powspctrm = TFR.powspctrm(index{:});
% and also correct time
TFR.time = TFR.time(index{timedim});

