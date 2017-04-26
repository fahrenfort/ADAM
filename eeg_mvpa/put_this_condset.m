function condSet = put_this_condset(condSet,thisCondSet,cSet)
% function condSet = put_this_condset(condSet,thisCondSet,cSet)
% puts new labels in condition set, used internally
% The condition set has the form condSet{class_nr}(set_nr,cond_label_nr)
% in which class_nr refers to the stimulus class, set_nr refers to the
% dataset (1 = train, 2 = test) and cond_label_nr refers to the original
% labels in the eeg file belonging to that stimulus class, so that 
% condSet{1} = [ 1, 2, 3; 99, 100, 101 ];
% condSet{2} = [ 4, 5, 6; 102, 103, 104];
% would mean that stimulus class 1 contains labels 1, 2 and 3 in the train
% set and 99, 100 and 101 in the test set, whereas stimulus class 2
% contains labels 4, 5 and 6 in the train set and 102, 103 and 104 in the
% test set. condSet can contain NaN's, but these are cut out whenever
% possible.
% if thisCondSet{1} = [1,2]; and thisCondSet{2} = [3,4];
% put_this_condset(condSet,thisCondSet,2) would replace [99, 100, 101] with
% [1,2] and replace [102, 103, 104] with [3,4], filling up the gaps with
% NaNs
indx = [2 1]; % get the other set
oSet = indx(cSet);
for c = 1:numel(condSet)
    condSet{c}(:,isnan(condSet{c}(oSet,:))) = []; % cut out NaNs right away
    orS = size(condSet{c},2);
    thisS = numel(thisCondSet{c});
    if thisS > orS % fill up with NaNs
        condSet{c}(cSet,orS+1:thisS) = NaN;
        condSet{c}(oSet,orS+1:thisS) = NaN;
    end
    if orS > thisS
        thisCondSet{c}(thisS+1:orS) = NaN;
    end
    condSet{c}(cSet,:) = thisCondSet{c};
end