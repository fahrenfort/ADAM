function [BDM, FEM] = EEG_backward_and_forward_matrix(allTestData,allTrainData,labelsOfTestSet,labelsOfTrainSet,msettings)
% function [BDM, FEM] = EEG_backward_and_forward_matrix(allTestData,allTrainData,labelsOfTestSet,labelsOfTrainSet,msettings)
% allTestData contains an array with conditions for testing, with
% dimensions (electrode x time x trial)
% Similarly, allTrainData contains an array with conditions for testing,
% with dimensions (electrode x time x trial)
% labelsOfTestSet and labelsOfTrainSet are 1-D arrays containing
% the labels corresponding to the conditions in the 3rd dimension of the
% train and test data.
% msettings is a struct containing some fields such as crossclass.
% if crossclass is 1, a crossclassification over time is executed
% otherwise, training and testing is only performed along the diagonal
% Method is the classification method: 'linear' (default), 'diagLinear' etc
% Type: help classify to see the options
% LabelsOverTime(t1,t2,actualLabel,classifierLabel) contains a
% frequency table for each time-time point, containing the frequency with
% which a classifierLabel was assigned to each stimulus type (i.e.
% actualLabel).
% Percentage correct can flexibly be computed from that table, e.g. by
% simply taking the number of times a label was assigned correctly
% WeightsOverTime(t,:) contains the electrode weights over time
% covPatternsOverTime(t,:) contains the weights multiplied by the covariance
% matrix to produce patterns that can be used for source localization
% Function can do both backward decoding (BDM) and forward encoding
% modeling (FEM). 

% J.J.Fahrenfort, VU 2015, 2016

% initialize input
doFEM = true;
doBDM  = true;
labelsonly = false;
method = 'linear';
crossclass = 0;
basis_sigma = 1;
v2struct(msettings);

% specify sigma if not already specified
nCond = numel(unique(labelsOfTrainSet));

% initialize output
FEM = [];
BDM = [];

% check input
if nCond ~= numel(unique(labelsOfTestSet)) && ~labelsonly
    error('number of conditions in train and test not equal, this only works if you specify method = ''labelsonly''');
end

% in case of FEM, create basis set
if doFEM
    if basis_sigma == 0
        disp('using a delta function (box car) as basis set');
        basis_set = eye(nCond,nCond); % channel * condition: the most basic basis set, any other set would impose a tuning curve on the data
    else
        disp(['using a gaussian function with a sigma of ' num2str(basis_sigma) ' as basis set']);
        set = mygausswin(1:nCond,basis_sigma,round(nCond/2))'; % use gaussian tuning curve, center on middle
        % set = set/sum(set); % normalize so that the gaussian sums to 1
        set = set./max(set); % normalize the max to 1
        basis_set = nan(nCond,nCond);
        for c=1:nCond
            basis_set(:,c) = circshift(shiftdim(set),-floor(nCond/2)+c); % make predictors for all channels
        end
    end
    FEM.basis_set = basis_set;
    % and pre-allocate some space for FEM
    C1 = nan(size(allTrainData,3),nCond); % trials * channel response
    FEM.WeightsOverTime = nan(size(allTrainData,2),size(allTrainData,1),nCond); % time * electrode * channels
    FEM.C2_percondition = nan(size(allTrainData,2),size(allTestData,2),nCond,nCond); % time * time * conditions * channel responses
    FEM.C2_average = nan(size(allTrainData,2),size(allTestData,2),nCond); % time * time * channel response
    for c=1:nCond
        indx = labelsOfTrainSet == c;
        C1(indx,:) = repmat(basis_set(:,c)',[sum(indx) 1]); % put basis set: C1 = trials * channel_response
    end
end

% in case of BDM, pre-allocate some space
if doBDM
    if labelsonly
        disp('agnostic about category labels in the test data, only returning classifier labels for all time points in each trial, no accuracy scores or aggregated response matrix data are computed');
        BDM.LabelsOverTime = zeros(size(allTestData,3),size(allTestData,2)); % trial * time
    else
        BDM.LabelsOverTime = zeros(size(allTestData,2),size(allTrainData,2),nCond,nCond); % time * time * cond * cond confusion matrix
    end
    BDM.WeightsOverTime = nan(size(allTrainData,2),size(allTrainData,1));
    BDM.covPatternsOverTime = nan(size(BDM.WeightsOverTime));
    BDM.corPatternsOverTime = nan(size(BDM.WeightsOverTime));
end

if ~verLessThan('matlab','8.3')
    disp('Running the faster classification algorithm (from Matlab R2014a onwards).');
end

% run through time * time
t2start = 1;
t2stop = size(allTestData,2);
% train loop
for t=1:size(allTrainData,2)
    fprintf(1,'.');
    if ~mod(t,80)
        fprintf(1,'\n');
    end
    dataTrain= squeeze(allTrainData(:,t,:))'; % transposed so we get trial x electrode instead of electrode * trial
    % C1 = trials * channel responses, B1 = trials * electrodes (dataTrain == B1)
    
    % get FEM weights magic: solve C1\B1 
    if doFEM 
        FEM_weights = (C1\dataTrain)'; % output transposed so FEM_weights = electrodes * channel responses
        FEM.WeightsOverTime(t,:,:) = FEM_weights; % time * electrodes * channel_responses
    end
    % get BDM weights, use fitcdiscr from R2014a onwards
    if doBDM && ~verLessThan('matlab','8.3')
        class_obj = compact(fitcdiscr(dataTrain,labelsOfTrainSet,'DiscrimType',method));
    end
    % test loop (or diagonal)
    if ~crossclass
        t2start = t;
        t2stop = t;
    end
    for t2=t2start:t2stop
        % get data
        dataTest=squeeze(allTestData(:,t2,:))'; % transpose so we get trial x electrode instead of electrode * trial
        % BDM: do test label matrix using backward model
        if doBDM
            % this is much more efficient when not "re-training" for every new test set
            % use fitcdiscr from R2014a onwards:
            if ~verLessThan('matlab','8.3')
                assignedLabels = predict(class_obj,dataTest);
                coeffs = class_obj.Coeffs;
            else
                [assignedLabels,~,~,~,coeffs] = classify(dataTest,dataTrain,labelsOfTrainSet,method);
            end
            if ~exist('coeffs','var') || isempty(coeffs)
                error('Classify is not returning coefficients, maybe you are not using the matlab native classify function. Check whether you have the biosig plugin installed in eeglab: this toolbox contains a classify function that you should not be using -> remove it or remove the path that leads to it.');
            end
            if labelsonly % just returning labels, only the diagonal, no cross corr (or we will run out of memory)
                BDM.LabelsOverTime(:,t) = assignedLabels;
            else % or the response matrix
                labelMatrix = zeros(nCond,nCond);
                for cActualLabel=1:nCond
                    for cClassifierLabel=1:nCond
                        labelMatrix(cActualLabel,cClassifierLabel) = numel(find(labelsOfTestSet==cActualLabel & assignedLabels==cClassifierLabel));
                    end
                end
                BDM.LabelsOverTime(t2,t,:,:) = labelMatrix; % time * time * cond * cond confusion matrix
            end
        end
        
        % FEM: do forward encoding model
        if doFEM
            % get hypothetical channel responses, main magic, solve w\B2 -> C2
            % FEM_weights = electrodes * channel responses, B2 = trial * electrodes (dataTest == B2)
            C2 = FEM_weights\(dataTest'); % output: channel response * trials (B2 is transposed to get electrodes * trial)
            % compute average channel response per condition
            C2_percondition = nan(nCond,nCond);
            C2_aligned = nan(nCond,nCond);
            for c = 1:nCond
                C2_percondition(c,:) = mean(C2(:,labelsOfTestSet==c),2); % cond * channel_response
                C2_aligned(c,:) = circshift(shiftdim((C2_percondition(c,:))),floor(nCond/2)-c); % align -> shift back to old location
            end
            % compute average C2 and retain outcomes
            FEM.C2_average(t2,t,:) = mean(C2_aligned,1); % average channel responses, time * time * channel_response
            FEM.C2_percondition(t2,t,:,:) = C2_percondition; % time * time * cond * channel_response
            % compute labels
            assignedLabels = nan(size(C2,2),1);
            for cTrial = 1:size(C2,2)
                if basis_sigma == 0
                    [~, assignedLabels(cTrial)] = max(C2(:,cTrial)); % just pick the peak
                else
                    [~, assignedLabels(cTrial)] = max(corr(C2(:,cTrial),basis_set)); % pick the condition for which the generated channel reponse correlates best
                end
            end
            % get response matrix
            FEM_labelMatrix = zeros(nCond,nCond);
            for cActualLabel=1:nCond
                for cClassifierLabel=1:nCond
                    FEM_labelMatrix(cActualLabel,cClassifierLabel) = numel(find(labelsOfTestSet==cActualLabel & assignedLabels==cClassifierLabel));
                end
            end
            FEM.LabelsOverTime(t2,t,:,:) = FEM_labelMatrix; % time * time * cond * cond confusion matrix
        end
        
    end
    
    % weights for backward model (note that weights depend on training set only)
    % also compute the pattern matrix from the weights by multiplying
    % the covariance matrix with the classifier weights (see Haufe et al, Neuroimage, 2014)
    % dataTrain should be defined as trial x electrode
    
    % suggested modification Joram 07/14/17
    % in case of three (or more) categories, the weights used to only give
    % the separability of the first two categories;
    % below I loop over the possible combinations, and average the weights
    % over these combinations; loop should also work for two categories,
    % where it just goes over the loop once to get coeffs(1,2).Linear like
    % in the previous code.
    if doBDM
        if ~verLessThan('matlab','8.3') % again, slightly different for more recent version of matlab
            k=1;
            for r=1:length(coeffs)-1
                for c = r+1:length(coeffs)
                    weights(:,k) = shiftdim(coeffs(r,c).Linear);
                    k=k+1;
                end
            end
        else
            k=1;
            for r=1:length(coeffs)-1
                for c = r+1:length(coeffs)
                    weights(:,k) = shiftdim(coeffs(r,c).linear);
                    k=k+1;
                end
            end
        end
        BDM.WeightsOverTime(t,:) = mean(weights,2); 
        BDM.covPatternsOverTime(t,:) = cov(dataTrain)*mean(weights,2);
        BDM.corPatternsOverTime(t,:) = corr(dataTrain)*mean(weights,2);
    end
    
end
fprintf(1,'\n');
