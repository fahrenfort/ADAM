function measure = class_accuracy_from_matrix(labelMatrix,method,crossclass)
% function measure = class_accuracy_from_matrix(labelMatrix,method,crossclass)
% Generates an accuracy measure from a confusion matrix LabelMatrix
% (actual,classified). Labelmatrix can also be a crosscorrelation matrix of
% LabelMatrix over time: labelMatrix(time1, time2, actual,classified)
% contains frequency counts for each combination of the category label
% (actual) and the assigned label by the response system (classified) for
% each time_train time_test combination.
% Available methods are:
% 'accuracy', which computes the average number of correct
% classifications over all category labels (default)  
% 'hr-far', which computes the accuracy as the hit rate (number of
% correctly classified items over all items of category present) minus the
% false alarm rate (number of times that category was falsely identified
% when the category was absent), assuming the first element contains counts
% of the true signal and the second element contains counts of the noise
% signal.
% 'dprime' has the same requirements as 'hr-far' but computes d'.
% 'hr' just gives hit rate
% 'far' just gives false alarm rate
% 'mr' just gives miss rate
% 'cr' just gives correct rejection rate
% LabelMatrix can also be a matrix over time: labelMatrixOverTime(t,t2,actual,classified)
% in which case the function computes the output for all t, t2 time
% combinations.
%
% By J.J.Fahrenfort, VU, 2015

if nargin < 3
    crossclass = true;
end
if nargin < 2
    method = 'accuracy';
end
if ismatrix(labelMatrix) % does not contain time
    measure = class_accuracy_from_matrix_only(labelMatrix,method);
elseif ndims(labelMatrix) == 4 % wrapper in case labelMatrix also contains time as the first two dimensions
    measure = zeros(size(labelMatrix,1), size(labelMatrix,2));
    t2start = 1;
    t2stop = size(labelMatrix,2);
    for t=1:size(labelMatrix,1) % test loop
        if ~crossclass % train loop (or diagonal)
            t2start = t;
            t2stop = t;
        end
        for t2=t2start:t2stop
            measure(t,t2) = class_accuracy_from_matrix_only(squeeze(labelMatrix(t,t2,:,:)),method);
        end
    end
else
    error('the labelMatrix array does not contain a sensible number of dimensions');    
end

function measure = class_accuracy_from_matrix_only(labelMatrix,method)

% make sure we get a symmetrical matrix
if size(labelMatrix,1) ~= size(labelMatrix,2) || ~ismatrix(labelMatrix)
    error('the response count matrix is not symmetrical: it should have actual labels as rows and response labels as columns');
end
nCond = size(labelMatrix,1);

% compute accuracy measure
if strcmpi(method,'accuracy')
    percCorrect = zeros(nCond,1);
    for cActualLabel = 1:nCond
        nActualLabel = sum(labelMatrix(cActualLabel,:));
        percCorrect(cActualLabel) = labelMatrix(cActualLabel,cActualLabel)/nActualLabel;
    end
    measure = mean(percCorrect);
elseif any(strcmpi(method,{'hr-far','dprime','hr','far','mr','cr'}))
    if nCond ~= 2
        error('cannot compute hr-far because the number of stimulus categories is unequal to 2');
    end
    % disp('important: assuming the hr/far should be computed on stimulus present label = 1, stimulus absent label = 2');
    stimPresent = 1;
    stimAbsent = 2;
    if strcmpi(method,'dprime')
        % correct for infinity using Hautus MJ (1995) Corrections for
        % extreme proportions and their biasing effects on estimated values
        % of d?. Behavior Research Methods, Instruments, & Computers 27(1):46?51
        labelMatrix = labelMatrix + .5;
    end
    nStimPresent = sum(labelMatrix(stimPresent,:));
    hr = labelMatrix(stimPresent,stimPresent)/nStimPresent;
    nStimAbsent = sum(labelMatrix(stimAbsent,:));
    far = labelMatrix(stimAbsent,stimPresent)/nStimAbsent;
    if strcmpi(method,'hr-far')
        measure = hr - far;
    elseif strcmpi(method,'dprime')
        measure = dprime(hr,far);
    elseif strcmpi(method,'hr')
        measure = hr;
    elseif strcmpi(method,'far')
        measure = far;
    elseif strcmpi(method,'mr')
        measure = 1-hr;
    elseif strcmpi(method,'cr')
        measure = 1-far;
    end
else
    error('no existing method was specified');
end