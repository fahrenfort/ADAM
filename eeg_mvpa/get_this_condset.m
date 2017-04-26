function thisCondSet = get_this_condset(condSet,cSet)
% function thisCondSet = get_this_condset(condSet,cSet)
% obtains labels from the condition set, used internally
% The condition set has the form condSet{class_nr}(set_nr,cond_label_nr)
% in which class_nr refers to the stimulus class, set_nr refers to the
% dataset (1 = train, 2 = test) and cond_label_nr refers to the original
% labels in the eeg file belonging to that stimulus class, so that 
% condSet{1} = [ 1, 2, 3; 99, 100, 101 ];
% condSet{2} = [ 4, 5, 6; 102, 103, 104];
% would mean that stimulus class 1 contains labels 1, 2 and 3 in the train
% set and 99, 100 and 101 in the test set, whereas stimulus class 2
% contains labels 4, 5 and 6 in the train set and 102, 103 and 104 in the
% test set. condSet can contain NaN's, but these are removed when using
% this function to obtain labels.
% get_this_condset(condSet,2) would return the test labels:
% thisCondset{1} = [ 99, 100, 101 ];
% thisCondset{2} = [ 102, 103, 104];
thisCondSet{numel(condSet)} = [];
for c = 1:numel(condSet)
    thisCondSet{c} = condSet{c}(cSet,~isnan(condSet{c}(cSet,:)));
end