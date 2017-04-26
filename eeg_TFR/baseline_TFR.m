function TFR = baseline_TFR(TFR, condSet, tf_baseline, bl_method, subtr_method)
% function TFR = baseline_TFR(TFR, condSet, tf_baseline, bl_method, subtr_method)
%
% Performs baselining per condition on time frequency data (this is
% particularly important if you expect baseline differences per condition,
% but are not interested in those pre-onset baseline differences)
% asumes the data is in format 'rpt_chan_freq_time'
%
% Johannes Fahrenfort, VU 2016

if ~strcmp(TFR.dimord,'rpt_chan_freq_time')
    error('data in wrong format, rpt_chan_time required');
end

% per original label, or per new class label?
if strcmpi(subtr_method,'subtr_indiv')
    condSet = num2cell(unique([condSet{:}]));
    disp('Applying TFR baseline correction on every original condition label');
else
    disp('Applying TFR baseline correction on every stimulus class');
end

% do this for each condition (can easily be changed if so desired)
for c=1:numel(condSet)
    data = TFR.powspctrm(ismember(TFR.trialinfo,condSet{c}),:,:,:);
    bl_timebool = TFR.time>=tf_baseline(1) & TFR.time<=tf_baseline(2);
    meandata = repmat(nanmean(data(:,:,:,bl_timebool), 4), [1 1 1, size(data, 4)]);
    
    % do actual baseline operation
    if (strcmpi(bl_method, 'absolute'))
        data = data - meandata;
    elseif (strcmpi(bl_method, 'relative'))
        data = data ./ meandata;
    elseif (strcmpi(bl_method, 'relchange'))
        data = (data - meandata) ./ meandata;
    elseif (strcmpi(bl_method, 'vssum'))
        data = (data - meandata) ./ (data + meandata);
    elseif (strcmpi(bl_method, 'db'))
        data = 10*log10(data ./ meandata);
    end
    
    % inject back into data structure
    TFR.powspctrm(ismember(TFR.trialinfo,condSet{c}),:,:,:) = data;
    
end