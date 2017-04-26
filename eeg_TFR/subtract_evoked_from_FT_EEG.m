function FT_EEG = subtract_evoked_from_FT_EEG(FT_EEG,FT_EVOKED_EEG,field)
% Takes FT_EEG and FT_EVOKED_EEG (bins) as input, subtracting the
% evoked responses from the corresponding single trials. field specifies
% the field on which to execute the subtraction (default: 'trial'). First
% use  compute_bins_on_FT_EEG or to compute_erp_on_FT_EEG compute the
% evoked set. Then run this function to compute the induced set. After
% this, compute TFR using compute_TFR_EEG, and finally run ft_freqbaseline
% before doing MVPA or other. 
% IMPORTANT: MAKE SURE YOU BASELINE FT_EEG BEFORE RUNNING THIS FUNCTION
%
% J.J.Fahrenfort 2014, VU

% start with a whole bunch of sanity checks
if nargin < 3
    field = 'trial';
end
if ~isfield(FT_EVOKED_EEG,'oldindex')
    error('can not run subtraction procedure on these datasets, no information on how to match up trials.');
end
if ~isfield(FT_EEG,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on FT_EEG to resolve this.');
else
    dims = regexp(FT_EEG.dimord, '_', 'split');
    chandim = find(strcmp(dims,'chan'));
    timedim = find(strcmp(dims,'time'));
    trialdim = find(strcmp(dims,'rpt'));
end
if ~isfield(FT_EVOKED_EEG,'dimord')
    error('the input dataset is not in standard fieldtrip format, missing dimord.');
else
    if ~strcmp(FT_EVOKED_EEG.dimord,'rpt_chan_time');
        error('incorrect dimensions in FT_EVOKED_EEG: should have trial, channel and time (in that order) but not freq');
    end
end
if isempty(chandim | timedim | trialdim) || numel(dims) ~= 3
    error('incorrect dimensions of FT_EEG: should have trial, channel and time, but not frequency');
end
if ~isfield(FT_EEG,field) || ~isfield(FT_EVOKED_EEG,field)
    error([field ' does not seem to be contained in FT_EEG and/or FT_EVOKED_EEG']);
end
% make sure that electrode labels match up
chanindexFT = zeros(1,numel(FT_EEG.label));
chanindexEVK = zeros(1,numel(FT_EVOKED_EEG.label));
for cChan = 1:numel(FT_EEG.label)
    chanindexFT(cChan) = cChan;
    chanindexEVK(cChan) = find(strcmp(FT_EEG.label{cChan},FT_EVOKED_EEG.label));
end
% extract relevant data, in the correct dimorder
trial = permute(FT_EEG.(field),[trialdim chandim timedim]); 
evoked = FT_EVOKED_EEG.(field); % should already be 'rpt_chan_time'
% subtract the evoked responses from the single trials
for cBin=1:size(evoked,1)
    index = squeeze(FT_EVOKED_EEG.oldindex{cBin});
    trial(index,chanindexFT,:) = trial(index,chanindexFT,:) - repmat(evoked(cBin,chanindexEVK,:),[numel(index) 1 1]);
end
% remove trials without an evoked response and insert new data
index2remove = setdiff(1:size(trial,1),unique(vertcat(FT_EVOKED_EEG.oldindex{:})));
disp(['removing ' num2str(numel(index2remove)) ' trials because they either do not belong to a complete subset or were not epoched at all']);
trial(index2remove,:,:) = [];
FT_EEG.(field) = trial;
FT_EEG.trialinfo(index2remove) = []; % just keep the old trial info
FT_EEG.dimord = 'rpt_chan_time';
if isfield(FT_EEG,'sampleinfo');
    FT_EEG.sampleinfo(index2remove,:) = [];
end



