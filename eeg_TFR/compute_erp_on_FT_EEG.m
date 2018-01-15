function FT_EEG = compute_erp_on_FT_EEG(FT_EEG,condSet,field,method)
% function FT_EEG = compute_erp_on_FT_EEG(FT_EEG,condSet,field,method)
% Compute erp from fieldtrip standard format, necessary to compute
% induced single trial TFR data. The data format can only contain
% time, channel and trials. condSet can be set up as follows (example):
% condSet{1} = [ 1, 2, 3];
% condSet{2} = [ 4, 5, 6];
% averages single instances from condition 1,2 and 3 into a 'new' erp
% named condition 1, and averages single instances from condition 4, 5 and
% 6 into an erp of 'new' condition 2.
% If method is 'indiv' computes an erp seperately for every condition in condset 
% IMPORTANT: MAKE SURE YOU BASELINE FT_EEG BEFORE RUNNING THIS FUNCTION OR
% YOU MIGHT BE WASHING OUT YOUR ERPs
% J.J.Fahrenfort, VU 2014
if nargin < 4
    method = 'bin';   
end
if nargin < 3
    field = 'trial';
end
if nargin < 2
    error('missing condset when calling compute_erp_on_FT_EEG');
end
if size(condSet{1},1) > 1
    error('Pass only the condSet that belongs to this trial set to compute_erp_on_FT_EEG. Use get_this_condset.m to select the correct one');
end
if ~isfield(FT_EEG,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on FT_EEG to resolve this.');
else
    dims = regexp(FT_EEG.dimord, '_', 'split');
    chandim = find(strcmp(dims,'chan'));
    timedim = find(strcmp(dims,'time'));
    trialdim = find(strcmp(dims,'rpt'));
end

if isempty(chandim | timedim | trialdim) || numel(dims) > 3
    error('incorrect dimensions: should have time, channel and trial, cannot compute ERP on FT_EEG');
end
if ~isfield(FT_EEG,field)
    error([field ' should contain the data, but ' field ' is not a field of FT_EEG']);
end

% get trial info
trialinfo = FT_EEG.trialinfo;
trial = permute(FT_EEG.(field),[trialdim chandim timedim]);
if strcmp(method,'indiv')
    condSet = num2cell(unique([condSet{:}]));
end
for cCondSet = 1:numel(condSet)
    thisCondSet = [condSet{cCondSet}];
    if ischar(thisCondSet)
        thisCondSet = string2double(thisCondSet);
    end
    trialindex = find(ismember(trialinfo,thisCondSet));
    oldindex{cCondSet} = trialindex;
    finaltrialinfo(cCondSet,1) = cCondSet;
    finaltrial(cCondSet,:,:) = squeeze(mean(trial(trialindex,:,:),1));
end
FT_EEG.dimord = 'rpt_chan_time'; % should be 'rpt_chan_time' for consistency with conversion to raw format in ft_checkdata
FT_EEG.(field) = finaltrial;
FT_EEG.trialinfo = finaltrialinfo;
if isfield(FT_EEG,'oldindex')
    for c=1:numel(oldindex)
        oldindex{c} = {FT_EEG.oldindex{oldindex{c}}};
    end
end
FT_EEG.oldindex = oldindex;
if isfield(FT_EEG,'sampleinfo');
    FT_EEG = rmfield(FT_EEG,'sampleinfo');
end
if isfield(FT_EEG,'cfg');
    FT_EEG = rmfield(FT_EEG,'cfg');
end