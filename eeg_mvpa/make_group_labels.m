function grplabels = make_group_labels(trialinfo,condSet)
% creates group labels from condset and trial info
% J.J.Fahrenfort, VU, 2016
grplabels = zeros(numel(trialinfo),1);
for cCondSet=1:numel(condSet)
    grplabels(ismember(trialinfo,condSet{cCondSet})) = cCondSet;
end