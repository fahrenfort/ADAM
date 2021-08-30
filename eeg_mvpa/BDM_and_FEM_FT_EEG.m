function [BDM, FEM] = BDM_and_FEM_FT_EEG(train_FT_EEG,test_FT_EEG,train_condSet,test_condSet,msettings)
% BDM_and_FEM_FT_EEG computes the backward decoding (BDM) and forward encoding model (FEM) based on
% a train-test procedure.
%
% INPUTS:
%
% train_FT_EEG      train_FT_EEG.trial contains an array with conditions for testing, with
%                   dimensions (channel * time * trial)
% test_FT_EEG       test_FT_EEG.trial contains an array with conditions for testing, with dimensions
%                   (channel * time * trial)
% train_condSet     class definitions containing the train labels corresponding to the classes
%                   train_FT_EEG.trialinfo
% test_condSet      class definitions containing the test labels corresponding to the classes
%                   test_FT_EEG.trialinfo
% msettings         a struct containing some fields such as crossclass. if crossclass is 1, a
%                   crossclassification over time is executed otherwise, training and testing is
%                   only performed along the diagonal. Method is the classification method: 'linear'
%                   (default), 'diagLinear' etc Type: help classify to see the options
%
% OUTPUTS:
%
% BDM and FEM are backward decoding and forward encodming models respectively, each containing
% fields such as LabelsOverTime(t1,t2,actualLabel,classifierLabel) contains a frequency table for
% each time-time point, containing the frequency with which a classifierLabel was assigned
% to each stimulus type (i.e. actualLabel). Percentage correct can flexibly be computed from that
% table, e.g. by simply taking the number of times a label was assigned correctly
% WeightsOverTime(t,:) contains the electrode weights over time covPatternsOverTime(t,:) contains
% the weights multiplied by the covariance matrix to produce patterns that can be used for source
% localization Function can do both backward decoding (BDM) and forward encoding modeling (FEM).
%
% Internal function of the ADAM toolbox by J.J.Fahrenfort, VU 2015, 2016, 2018
%
% See also: ADAM_MVPA_FIRSTLEVEL, CLASSIFY_RAW_EEGLAB_DATA, CLASSIFY_TFR_FROM_EEGLAB_DATA

% default parameters
doFEM = true;
doBDM  = true;
compute_performance = true;
save_confidence = false;
method = 'linear';
crossclass = 0;
basis_sigma = 1;
v2struct(msettings);

% fix dimord
train_FT_EEG = fix_dimord(train_FT_EEG,'chan_time_rpt'); % should be channel * time * trial
test_FT_EEG = fix_dimord(test_FT_EEG,'chan_time_rpt'); % should be channel * time * trial

% initialize input
nClasses = numel(train_condSet);
allTrainData = train_FT_EEG.trial;
labelsOfTrainSet = make_group_labels(train_FT_EEG.trialinfo, train_condSet);
allTestData = test_FT_EEG.trial;
if compute_performance
    % check input
    if nClasses ~= numel(test_condSet)
        error('number of classes in train and test not equal, this only works if you do not compute performance but only classifier confidence');
    end
    % make class labels
    labelsOfTestSet = make_group_labels(test_FT_EEG.trialinfo, test_condSet);
else
    labelsOfTestSet =  NaN(size(test_FT_EEG.trialinfo)); % when no meaningful class definition exists in the test set
end

% initialize output
FEM = [];
BDM = [];

% in case of FEM, create basis set
if doFEM
    if basis_sigma == 0
        disp('using a delta function (box car) as basis set');
        basis_set = eye(nClasses,nClasses); % channel * condition: the most basic basis set, any other set would impose a tuning curve on the data
    else
        disp(['using a gaussian function with a sigma of ' num2str(basis_sigma) ' as basis set']);
        set = mygausswin(1:nClasses,basis_sigma,round(nClasses/2))'; % use gaussian tuning curve, center on middle
        % set = set/sum(set); % normalize so that the gaussian sums to 1
        set = set./max(set); % normalize the max to 1
        basis_set = nan(nClasses,nClasses);
        for c=1:nClasses
            basis_set(:,c) = circshift(shiftdim(set),-floor(nClasses/2)+c); % make predictors for all channels
        end
    end
    FEM.basis_set = basis_set;
    % and pre-allocate some space for FEM
    C1 = nan(size(allTrainData,3),nClasses); % trials * channel response
    FEM.WeightsOverTime = nan(size(allTrainData,2),size(allTrainData,1),nClasses); % time * electrode * channels
    FEM.C2_percondition = nan(size(allTrainData,2),size(allTestData,2),nClasses,nClasses); % time * time * conditions * channel responses
    FEM.C2_average = nan(size(allTrainData,2),size(allTestData,2),nClasses); % time * time * channel response
    for c=1:nClasses
        indx = labelsOfTrainSet == c;
        C1(indx,:) = repmat(basis_set(:,c)',[sum(indx) 1]); % put basis set: C1 = trials * channel_response
    end
end

% in case of BDM, pre-allocate some space
if doBDM
    BDM.WeightsOverTime = NaN(size(allTrainData,2),size(allTrainData,1));
    BDM.covPatternsOverTime = NaN(size(BDM.WeightsOverTime));
    BDM.corPatternsOverTime = NaN(size(BDM.WeightsOverTime));
    if save_confidence
        BDM.TrueClassLabels = labelsOfTestSet;
    end
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
        if strcmpi(method,'svm')
            class_obj = compact(fitcecoc(dataTrain,labelsOfTrainSet));
            coeffs = NaN;
        else
            class_obj = compact(fitcdiscr(dataTrain,labelsOfTrainSet,'DiscrimType',method));
            coeffs = class_obj.Coeffs;
        end
    end
    % test loop (or diagonal)
    if ~crossclass
        t2start = t;
        t2stop = t;
    end
    for t2=t2start:t2stop
        
        % get data
        dataTest=squeeze(allTestData(:,t2,:))'; % transpose so we get trial x electrode instead of electrode * trial
        
        % BDM: do backward decoding model
        if doBDM
            % this is much more efficient when not "re-training" for every new test set
            % use fitcdiscr from R2014a onwards:
            if ~verLessThan('matlab','8.3')
                [~,scores,~]  = predict(class_obj,dataTest);
            else
                [~,~,scores,~,coeffs] = classify(dataTest,dataTrain,labelsOfTrainSet,method);
            end
            if ~exist('coeffs','var') || isempty(coeffs) && ~strcmpi(method,'svm')
                error('Classify is not returning coefficients, maybe you are not using the matlab native classify function. Check whether you have the biosig plugin installed in eeglab: this toolbox contains a classify function that you should not be using -> remove it or remove the path that leads to it.');
            end
            % compute performance
            if compute_performance
                BDM.ClassOverTime(t2,t) = compute_class_performance(scores,labelsOfTestSet,measuremethod);
            else
                BDM.ClassOverTime(t2,t) = mean(scores(:,1),1); % return the probability of assigning the label of the first class
            end
            % save classifier confidence scores
            if save_confidence && t == t2 % also returning confidence scores, Only the diagonal or we will run out of memory, but can be changed if necessary.
                BDM.ScoresOverTime(:,t,:) = scores;
            end
        end
        
        % FEM: do forward encoding model
        if doFEM
            % get hypothetical channel responses, main magic, solve w\B2 -> C2
            % FEM_weights = electrodes * channel responses, B2 = trial * electrodes (dataTest == B2)
            C2 = FEM_weights\(dataTest'); % output: channel response * trials (B2 is transposed to get electrodes * trial)
            % compute average channel response per condition
            C2_percondition = nan(nClasses,nClasses);
            C2_aligned = nan(nClasses,nClasses);
            for c = 1:nClasses
                C2_percondition(c,:) = mean(C2(:,labelsOfTestSet==c),2); % cond * channel_response
                C2_aligned(c,:) = circshift(shiftdim((C2_percondition(c,:))),floor(nClasses/2)-c); % align -> shift back to old location
            end
            % compute average C2 and retain outcomes
            FEM.C2_average(t2,t,:) = mean(C2_aligned,1); % average channel responses, time * time * channel_response
            FEM.C2_percondition(t2,t,:,:) = C2_percondition; % time * time * cond * channel_response
        end
        
    end % end t2 loop
    
    % Compute weights for backward model (note that weights here depend on training set only) also
    % compute the pattern matrix from the weights by multiplying the covariance matrix with the
    % classifier weights (see Haufe et al, Neuroimage, 2014) dataTrain should be defined as trial x
    % electrode
    
    % Suggested modification Joram 07/14/17 in case of three (or more) categories, the weights used
    % to only give the separability of the first two categories; below I loop over the possible
    % combinations, and average the weights over these combinations; loop should also work for two
    % categories, where it just goes over the loop once to get coeffs(1,2).Linear like in the
    % previous code.
    % Remark JJF: fine, we'll leave it this way for now, but note that a topomap visualization of a
    % difference between three classes/conditions is simply not that informative. In that case it is
    % better to run three two-class contrasts to visualize the differences.
    if doBDM
        if ~strcmpi(method,'svm')
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
        end
        
        % to compute activation patterns we need the covariance of the original (unwhitened) data
        if ~exist('matObj','var')
            matObj = matfile(fname);
        end
        % slightly different for frequency or raw data
        if exist('dim_params','var') % frequency data
            dim_params.index{dim_params.timedim} = t; % only read for the current time point
            realData = squeeze(matObj.powspctrm(dim_params.index{:}));
            if dim_params.trialdim > dim_params.chandim
                realData = realData'; % transpose to make realData have trial x chan
            end
        else % raw data
            realData = squeeze(matObj.trial(:,:,t)); % realData is trial x chan
        end
        % finally compute activation pattern according to Haufe
        if ~strcmpi(method,'svm')
            BDM.covPatternsOverTime(t,:) = cov(realData)*mean(weights,2);
            BDM.corPatternsOverTime(t,:) = corr(realData)*mean(weights,2);
        end
        clear realData;
    end
    
end % end t1 loop

fprintf(1,'\n');
