function class_perf = compute_class_performance(scores,true_labels,measuremethod)
% compute_class_performance is an internal function of the ADAM toolbox
% It computes classifier performance, using the test scores of the classifier and the true labels of
% the test set. Replaces class_accuracy_from_matrix which is obsolete from ADAM 1.0.5 onwards.
%
% Inputs:
%       scores              confidences scores of the classifier to the trials in the test set
%       true_labels         true labels of the trials in the test set
%       measuremethod       the method used to compute classifier performance
%
% available measure methods are:
%       'AUC' (default) computes Area Under the Curve 
%       'accuracy'      computes balanced accuracy (number of correct classifications per class,
%                       averaged over all classes)
%       'hr-far'        computes the accuracy as the hit rate (number of correctly classified
%                       items over all items of category present) minus the false alarm rate (number
%                       of times that category was falsely identified when the category was absent).
%       'dprime'        computes d'.
%       'hr'            just gives hit rate
%       'far'           just gives false alarm rate
%       'mr'            just gives miss rate 
%       'cr'            just gives correct rejection rate
%
% Output: classifier performance
%
% By J.J.Fahrenfort, UvA/VU 2018
% See also: ADAM_MVPA_FIRSTLEVEL, CLASSIFY_RAW_EEGLAB_DATA, CLASSIFY_TFR_FROM_EEGLAB_DATA

% make sure we get a symmetrical matrix
nClasses = size(scores,2);

% compute the confusion matrix (not necessary for AUC)
if ~strcmpi(measuremethod,'AUC')
    confusionMatrix = zeros(nClasses,nClasses);
    [~, assigned_labels] = max(scores,[],2);
    for cActualLabel=1:nClasses
        for cClassifierLabel=1:nClasses
            confusionMatrix(cActualLabel,cClassifierLabel) = numel(find(true_labels==cActualLabel & assigned_labels==cClassifierLabel));
        end
    end
end

% compute classfier performance measure
if strcmpi(measuremethod,'AUC')
    % [~,~,~,AUC] = perfcurve(true_labels,scores(:,1),1); % THIS IS SLOW!
    % solve this for multi-class, and do it fast
    pairs = nchoosek(1:nClasses,2); % pair-wise combinations of classes
    if size(pairs,1)>1
        % do this both ways in case of multi-class problems
        pairs(end+1:end+size(pairs,1),:) = pairs(:,[2 1]);
    end
    AUC = zeros(1,size(pairs,1));
    for c=1:size(pairs,1)
        % grab all scores of the pair, the first one is always poslabel
        ind2use = ismember(true_labels,pairs(c,:));     % grab only two classes
        boollabels = false(size(true_labels));          % set all labels to false
        boollabels(true_labels==pairs(c,1)) =  true;    % set only positive class to true
        labels2use = boollabels(ind2use);           % select pairwise labels
        scores2use = scores(ind2use,pairs(c,1));    % select pairwise scores
        AUC(c) = scoreAUC(labels2use,scores2use);   % compute AUC, much FASTER! :-)
    end
    class_perf = mean(AUC);
elseif strcmpi(measuremethod,'accuracy')
    percCorrect = zeros(nClasses,1);
    for cActualLabel = 1:nClasses
        nActualLabel = sum(confusionMatrix(cActualLabel,:));
        percCorrect(cActualLabel) = confusionMatrix(cActualLabel,cActualLabel)/nActualLabel;
    end
    class_perf = mean(percCorrect);
elseif any(strcmpi(measuremethod,{'dprime','hr-far','hr','far','mr','cr'}))
    if nClasses ~= 2
        error('cannot compute hr-far because the number of stimulus classes is unequal to 2');
    end
    % correct for infinity or minus infinity using:
    % Hautus MJ (1995) Corrections for extreme proportions and their biasing effects on estimated
    % values of d'. Behavior Research measuremethods, Instruments, & Computers 27(1):46-51
    if strcmpi(measuremethod,'dprime') && any(any(confusionMatrix == 0))
        confusionMatrix = confusionMatrix + .5;
    end
    % When computing hr/far/mr/cr, assuming that:
    % stimulus present label = 1
    % stimulus absent label = 2
    % (for the outcome on d' or hr-far it does not matter whether signal or noise goes first)
    stimPresent = 1;
    stimAbsent = 2;
    nStimPresent = sum(confusionMatrix(stimPresent,:));
    hr = confusionMatrix(stimPresent,stimPresent)/nStimPresent;
    nStimAbsent = sum(confusionMatrix(stimAbsent,:));
    far = confusionMatrix(stimAbsent,stimPresent)/nStimAbsent;
    if strcmpi(measuremethod,'dprime')
        class_perf = norminv(hr) - norminv(far);
    elseif strcmpi(measuremethod,'hr-far')
        class_perf = hr - far;
    elseif strcmpi(measuremethod,'hr')
        class_perf = hr;
    elseif strcmpi(measuremethod,'far')
        class_perf = far;
    elseif strcmpi(measuremethod,'mr')
        class_perf = 1-hr;
    elseif strcmpi(measuremethod,'cr')
        class_perf = 1-far;
    end
else
    error('ERROR: No existing measuremethod for computing the performance of the classifier was specified');
end


