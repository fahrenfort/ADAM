function [ FT_EEG, newCondSet ] = compute_bins_on_FT_EEG(FT_EEG,condSet,field,labeltype)
% Compute binned data from fieldtrip standard format, necessary to compute
% evoked and induced single trial TFR data, but also to increase signal to
% noise ratio for subsequent MVPA analyses. The data format can contain 
% either time, channel and trials, but may also include frequency (for TFR
% data). condSet can be set up as follows (example): 
% condSet{1} = [ 1, 2, 3];
% condSet{2} = [ 4, 5, 6];
% averages single instances from condition 1,2 and 3 into a 'new' instance
% of condition 1, and averages single instances from condition 4, 5 and 6
% into a single new instance of condition 2. 'Leftover' trials are
% discarded, e.g. if there are 8 instances of 1, 5 instance of condition 2
% and 7 instances of condition 3, then the 'new' condition 1 will contain 5
% trials, the leftover trials (3 'old' condition 1 and 2 'old' condition 3)
% are discarded. The new stimulus classes are either given new condition
% labels, counting from 1 to the number of bins (labeltype = 'newlabel',
% default), or are re-labeled using the first condition label from each bin
% (specify labeltype = 'original'), in which case the trials containing the
% averages from condSet{1} would get label 1, but the trials containing the
% average from condSet{2} in the above example would get label 4.
% You can also bin multiple trials from the same condition together: 
% condSet{1} = [ 1, 1, 1];
% condSet{2} = [ 2, 2, 2];
% averages 3 instances from condition 1 together into a single instance and
% averages 3 instances from condition 2 together into a single instance
% Again, leftover trials are discarded (if there were 10 trials of
% condition 1, the new set will contain 3 trial averaged instances of
% condition 1, and 1 'old' trial will be discarded.
%
% By J.J.Fahrenfort, VU 2014, 2015

if nargin < 4
    labeltype = 'newlabel';
end
if nargin < 3
    field = 'trial';
end
if ~isfield(FT_EEG,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on FT_EEG to resolve this.');
else
    dims = regexp(FT_EEG.dimord, '_', 'split');
    chandim = find(strcmp(dims,'chan'));
    timedim = find(strcmp(dims,'time'));
    trialdim = find(strcmp(dims,'rpt'));
    freqdim = find(strcmp(dims,'freq'));
end
if size(condSet{1},1) > 1
    error('you should only pass the condSet of the dataset you want to bin, use the get_this_condset.m function to restrict condSet before passing it to compute_bins_on_FT_EEG.m');
end

if isempty(chandim | timedim | trialdim) || numel(dims) > 4
    error('incorrect dimensions: should have time, channel and trial, cannot bin FT_EEG');
end
if isempty(freqdim) && numel(dims) == 4
    error('has four dimensions, but does not include frequency, cannot bin FT_EEG');
end
if ~isfield(FT_EEG,field)
    error([field ' should contain the data, but ' field ' is not a field of FT_EEG']);
end
if ~isempty(freqdim)
    nFreqs = size(FT_EEG.(field),freqdim);
    freqs = true;
else
    nFreqs = 1;
    freqs = false;
end

% get the data, has form channel, time, trial and (optional) frequency
if freqs
    alltrial = permute(FT_EEG.(field),[chandim timedim trialdim freqdim]);
else
    alltrial = permute(FT_EEG.(field),[chandim timedim trialdim]);
end
    
for cFreq=1:nFreqs
    % get the frequency
    if freqs
        trial = squeeze(alltrial(:,:,:,cFreq));
    else
        trial = alltrial;
    end
    % get trial info
    trialinfo = FT_EEG.trialinfo;
    finaltrial = [];
    finaltrialinfo = [];
    oldindex = {};
    for cCondSet = 1:numel(condSet)
        % for each condition set, compute an array with index numbers
        condset = condSet{cCondSet}(~isnan(condSet{cCondSet}));
        nConds = numel(condset);
        clear temptrial tempindex; % here temptrial and tempindex are emptied
        for cCond = 1:nConds
            cond = condset(cCond); % the particular condition we're working on
            index = find(trialinfo==cond); % indexes all trials of that condition
            % check whether we find anything
            if isempty(index)
                error(sprintf(['cannot find any trials of condition ' num2str(cond) ' in trial list: ' cond_string(unique(trialinfo)') '\nYou may have to lower nFolds for this analysis to work.']));
            end
            % determine whether I should bin multiple instances of this condition
            bincount = (condset==cond);
            nRepeats = sum(bincount);
            if nRepeats > 1
                % determine the how-manieth condition of this type
                order = zeros(1,nConds);
                order(bincount) = 1:nRepeats;
                index = index(order(cCond):nRepeats:end);
            end
            temptrial{cCond} = trial(:,:,index); % put all trials of this subcondition in a cell
            tempindex{cCond} = index; % nr of trials in (sub)condition
        end
        % determine the least number of trials across all cells, and cut all cells down to this length
        minEpochs = min(cellfun('size', temptrial, 3));
        temptrial = cellfun(@(temptrial) temptrial(:,:,1:minEpochs), temptrial, 'UniformOutput', false);
        tempindex = cellfun(@(tempindex) tempindex(1:minEpochs), tempindex, 'UniformOutput', false);
        dim = ndims(temptrial{1}); % should always be 3: channel, time, trialnumber!
        if dim ~= 3
            disp('WARNING: only one trial in a subcondition');
            dim = 3;
        end
        temptrial = cat(dim+1,temptrial{:}); % concatenate cells back to an array (along an extra dimension)
        tempindex = num2cell(cat(2,tempindex{:}),2); % transpose -> keep track of the index numbers to which these correspond in the 'old' trials
        if ndims(temptrial) == dim+1
            temptrial = mean(temptrial,dim+1); % compute average across the instances of subconditions in this condition
        else
            % disp('it seems there is only one subcondition in this condition, no averaging was performed');
        end
        % add all trials to the total
        finaltrial = cat(3,finaltrial,temptrial);
        finaltrialinfo = cat(1,finaltrialinfo,repmat(cCondSet,minEpochs,1)); % create new condition labels
        oldindex = cat(1,oldindex,tempindex); % and also keep track of the old index numbers
        newCondSet{cCondSet} = cCondSet;
    end
    % and put the data back
    if freqs
        alldata(:,:,:,cFreq) = finaltrial; 
    else
        alldata = finaltrial;
    end
end
if freqs 
    FT_EEG.dimord = 'rpt_chan_freq_time'; % was 'chan_time_rpt_freq' but should be 'rpt_chan_freq_time' for consistency with conversion to raw format in ft_checkdata
    FT_EEG.(field) = permute(alldata,[3 1 4 2]);
else
    FT_EEG.dimord = 'rpt_chan_time'; % was 'chan_time_rpt' but should be 'rpt_chan_time' for consistency with conversion to raw format in ft_checkdata
    FT_EEG.(field) = permute(alldata, [3 1 2]); 
end
if ~strcmp(labeltype,'newlabel')
    % disp('keep the first label of each bin as the new condition label for that binned stimulus class');
    temptrialinfo = finaltrialinfo;
    for c = 1:numel(condSet)
        finaltrialinfo(ismember(temptrialinfo,c)) = condSet{c}(1);
    end
end
% get the goodies
FT_EEG.trialinfo = finaltrialinfo;
FT_EEG.oldindex = oldindex;
if isfield(FT_EEG,'sampleinfo');
    FT_EEG = rmfield(FT_EEG,'sampleinfo');
end
if isfield(FT_EEG,'cfg');
    FT_EEG = rmfield(FT_EEG,'cfg');
end