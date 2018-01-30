function FT_EEG = balance_FT_EEG(FT_EEG,condSet,whitened,field)
% balance_FT_EEG achieves between-class balancing of EEG data contained in an FT_EEG struct by
% oversampling stimulus classes using the ADASYN algorithm
%
% First averages instances across trigger codes to create within-class balanced instances,
% subsequently uses those instances to generate new instances using the ADASYN algorithm. See ADASYN
% for more details.
%
% Internal function of the ADAM toolbox, by J.J.Fahrenfort, 2018
%
% See also: ADAM_MVPA_FIRSTLEVEL, ADASYN_TIME_SERIES

if nargin < 4
    field = 'trial';
end
if ~isfield(FT_EEG,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on FT_EEG to resolve this.');
end

% check whether field is there
if ~isfield(FT_EEG,field)
    error([field ' should contain the data, but ' field ' is not a field of FT_EEG']);
end

% dimord should be 'rpt_chan_time' in this function, could rewrite to make sure input format does not change to save memory
FT_EEG = fix_dimord(FT_EEG,'rpt_chan_time',field); 

% get trial info
trialinfo = FT_EEG.trialinfo;
nConds = numel(condSet);

% Determine the number of instances in each class
for cCondSet = 1:nConds
    thisCondSet = [condSet{cCondSet}];
    if ischar(thisCondSet)
        thisCondSet = string2double(thisCondSet);
    end
    nEachClass(cCondSet) = sum(ismember(trialinfo,thisCondSet)); % count instances
end

% Create binned averages for every class instance (retaining within class balancing) 
% the binned averages are used to create new instances
FT_EEG_BINNED = compute_bins_on_FT_EEG(FT_EEG,condSet,field,'original');
bintrialinfo = FT_EEG_BINNED.trialinfo;

% For all majority-minority combinations, oversample class instances in the minority class
% using ADASYNs SMOTE algorithm
data = FT_EEG_BINNED.(field);
[~, majIndex] = max(nEachClass);
minIndices = setdiff(1:nConds,majIndex);
for c=1:numel(minIndices)
    minIndex = minIndices(c);
    majClass = [condSet{majIndex}];
    minClass = [condSet{minIndex}];
    nReq = nEachClass(majIndex)-nEachClass(minIndex); % how many duplications are required for this class
    in_labels = [zeros(sum(ismember(bintrialinfo,majClass)),1); ones(sum(ismember(bintrialinfo,minClass)),1)];
    in_data = [data(ismember(bintrialinfo,majClass),:,:); data(ismember(bintrialinfo,minClass),:,:)];
    out_data = ADASYN_time_series(in_data, in_labels, nReq, [], [], whitened);
    if ~isempty(out_data)
        FT_EEG.trial = [FT_EEG.trial; out_data];
        FT_EEG.trialinfo = [FT_EEG.trialinfo; repmat(minClass',[nReq/numel(minClass) 1])];
    end
end

