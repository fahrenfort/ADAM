function [ trainIndx, testIndx, nFolds ] = make_folds(trialinfo1,trialinfo2,condSet,nFolds,labelsOnly)
% function to split up the dataset(s) into two sets, used for training and
% testing, this function is used internally by classify_RAW_eeglab_data and
% classify_TFR_from_eeglab_data
% trainIndx{nFolds,nCondSet} contains index numbers of trialinfo1 for training
% testIndx{nFolds,nCondSet} contains index numbers of trialinfo2 for testing
%
% J.J.Fahrenfort, VU, 2015, 2016

if nargin<6
    labelsOnly = false;
end

% How many condition sets are there?
nCondSet = numel(condSet);

% check if we can have this number of folds at all, otherwise reduce
if nFolds > 1
    min_count = [];
    for cCondSet=1:nCondSet
        thisCondSet = get_this_condset(condSet,1);
        trainBool = ismember(trialinfo1,thisCondSet{cCondSet});
        thisCondSet = get_this_condset(condSet,2);
        testBool = ismember(trialinfo2,thisCondSet{cCondSet});
        commonIndex = intersect(find(trainBool),find(testBool));
        uniqueTrain = setdiff(find(trainBool),commonIndex);
        uniqueTest = setdiff(find(testBool),commonIndex);
        if ~isempty(commonIndex)
            min_count = min([min_count numel(commonIndex)]);
        end
        if ~isempty(uniqueTrain)
            min_count = min([min_count numel(uniqueTrain)]);
        end
        if ~isempty(uniqueTest)
            min_count = min([min_count numel(uniqueTest)]);
        end
    end
    if nFolds > min_count
        nFolds = min_count;
        disp(['WARNING: you specified more folds than is feasible given the number of trials you have and your stimulus class definitions. Lowering nFolds to: ' num2str(nFolds)]);
    end
end

% Create datasets for testing and training
trainIndx{nFolds,nCondSet} = []; % training
testIndx{nFolds,nCondSet} = []; % testing
for cCondSet=1:nCondSet
    % get training trials
    thisCondSet = get_this_condset(condSet,1);
    trainBool = shiftdim(ismember(trialinfo1,thisCondSet{cCondSet}));
    if isempty(find(trainBool,1))
        error(['Cannot find trials of one or more of the conditions in ' cond_string(thisCondSet{cCondSet}) ' in the training list: ' cond_string(unique(trialinfo1)')]);
    end
    % get testing trials
    thisCondSet = get_this_condset(condSet,2);
    testBool = shiftdim(ismember(trialinfo2,thisCondSet{cCondSet}));
    if labelsOnly
        testBool = ones(size(trialinfo2)); % pick all trials
    end
    if ~any(testBool)
        error(['Cannot find trials of one or more of the conditions in ' cond_string(thisCondSet{cCondSet}) ' in the testing list: ' cond_string(unique(trialinfo2)')]);
    end
    % split up into folds for the leave-one-out procedure (if applicable)
    if nFolds > 1 
        if numel(trainBool) ~= numel(testBool)
            error('You have traveled a wormhole, this cannot happen: your training and testing datasets should have been the same size given that they come from the same file, but they are not! AAAAARRRRGGGH!');
        end
        commonIndex = intersect(find(trainBool),find(testBool));
        uniqueTrain = setdiff(find(trainBool),commonIndex);
        uniqueTest = setdiff(find(testBool),commonIndex);
        % create folds, separately for common and unique so that their
        % relative contribution to each fold remains constant
        if ~isempty(commonIndex)
            if numel(commonIndex)/nFolds < 1
                error(['You are trying to run using ' num2str(nFolds) ' folds, but there are not enough trials in each condition to be able to do this. Lower nFolds when running this function.']);
            end
            borders = round(0:numel(commonIndex)/nFolds:numel(commonIndex));
        end
        if ~isempty(uniqueTrain)
            if numel(uniqueTrain)/nFolds < 1
                error(['You are trying to run using ' num2str(nFolds) ' folds, but there are not enough trials in each condition to be able to do this. Lower nFolds when running this function.']);
            end
            borders1 = round(0:numel(uniqueTrain)/nFolds:numel(uniqueTrain));
        end
        if ~isempty(uniqueTest)
            if numel(uniqueTest)/nFolds < 1
                error(['You are trying to run using ' num2str(nFolds) ' folds, but there are not enough trials in each condition to be able to do this. Lower nFolds when running this function.']);
            end
            borders2 = round(0:numel(uniqueTest)/nFolds:numel(uniqueTest));
        end
        for cFld=1:nFolds
            testIndex = [];
            trainIndex = [];
            if ~isempty(commonIndex)
                testIndex = commonIndex(borders(cFld)+1:borders(cFld+1));
                trainIndex = setdiff(commonIndex, testIndex);
            end
            if ~isempty(uniqueTrain) % take subset of size (1-nFolds)/nFolds 
                trainSubset = setdiff(uniqueTrain,uniqueTrain(borders1(cFld)+1:borders1(cFld+1)));
                trainIndex = sort([trainIndex; trainSubset]);
            end
            if ~isempty(uniqueTest) % take subset of size 1/nFolds
                testSubset = uniqueTest(borders2(cFld)+1:borders2(cFld+1));
                testIndex = sort([testIndex; testSubset]);
            end
            trainIndx{cFld,cCondSet} = trainIndex;
            testIndx{cFld,cCondSet} = testIndex;
        end
    else % or keep it simple without folding
        trainIndx{1,cCondSet} = find(trainBool);
        testIndx{1,cCondSet} = find(testBool);
    end
end