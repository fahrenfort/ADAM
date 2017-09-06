function [out_featuresSyn, out_labelsSyn] = ADASYN(in_features, in_labels, in_beta, in_kDensity, in_kSMOTE, in_featuresAreNormalized)
%this function implements the ADASYN method as proposed in the following
%paper:
%
%[1]: H. He, Y. Bai, E.A. Garcia, and S. Li, "ADASYN: Adaptive Synthetic
%Sampling Approach for Imbalanced Learning", Proc. Int'l. J. Conf. Neural
%Networks, pp. 1322--1328, (2008).
%
%the implementation follows the notation and equation numbers given in
%section 3.1.4 of another paper:
%
%[2]: H. He and E.A. Garcia, "Learning from imbalanced data",
%Knowledge and Data Engineering, IEEE Transactions on 21, no. 9,
%pp. 1263--1284, (2009).
%
%
%the purpose of the ADASYN method is to improve class balance towards
%equally-sized classes for a given input dataset. this is achieved by
%synthetically creating new examples from the minority class via linear
%interpolation between existing minority class samples. this approach is
%known as the SMOTE method, cf. section 3.1.3 in [2]. ADASYN is an
%extension of SMOTE, creating more examples in the vicinity of the boundary
%between the two classes, than in the interior of the minority class.
%cf. the supplied script demo_ADASYN for an example of this.
%
%
% INPUTS:
%----------
%in_features:
%(N \times P) matrix of numerical features. each row is one example, each
%column is one feature, hence there are N examples with P features each.
%
%in_labels:
%boolean N-vector of labels, defining the classes to which the examples in
%in_features belong.
%
%in_beta [default: 1]:
%desired level of balance, where 0 means that the size of the minority
%class will not be changed, and 1 means that the minority class will be
%ADASYNed to have (approximately, due to rounding) the same size as the
%majority class. any value of in_beta between 0 and 1 provides a compromise
%between these two extremes.
%note that in_beta IS NOT the resulting class ratio, but a percentage of
%how much class balance is improved in comparison to the given class
%balance! 0 means nothing is improved in comparison to the given class
%balance and 1 means class sizes are perfectly equalized (except for small
%rounding-related deviations).
%
%in_kDensity [default: 5]:
%k for kNN used in ADASYN density estimation, i.e. in calculation of the
%\Gamma_i values in eq. (4) of reference [2]. this is the kNN call that
%regards examples from both classes.
%
%in_kSMOTE [default: 5]: 
%k for kNN used in subsequent SMOTE-style synthesis of new examples.
%this is the kNN call that regards only examples from the minority class.
%cf. eq. (1) in reference [2].
%
%in_featuresAreNormalized [default: true]:
%boolean indicating whether the features (i.e. the different columns) in
%in_features are already normalized to the same scale or not.
%by default normalized features are assumed as the input, i.e. the user is
%expected to apply a normalization method of choice before passing the data
%to the ADASYN function.
%the practical difference in the two values of in_featuresAreNormalized is
%the following:
%true:  Euclidean distance is used in all kNN calls. for reasonable
%       results, in_features should already be normalized when calling
%       ADASYN().
%false: standardized Euclidean distance (type "doc knnsearch" into MATLAB
%       console and look for 'seuclidean' for an explanation) is used in
%       all kNN calls. any normalization already present in in_features is
%       ignored, and instead unit variance normalization is used.
%       however, this does NOT modify the data. the normalization is only
%       applied within knnsearch.
%
%
% OUTPUTS:
%----------
%out_featuresSyn, out_labelsSyn:
%features and labels of ONLY the synthetically created examples.
%note that each entry of out_labelsSyn is the label of the minority class
%since only examples of the minority class are created.
%concatenating [in_features out_featuresSyn] and [in_labels out_labelsSyn]
%gives a new example set with the desired class balance.
%
%
%
%
%-------------------------------------------------------------------------
% Version: 1.0
% Date: 2015-04-17
% Author: Dominic Siedhoff
%-------------------------------------------------------------------------
%
%
%-------------------------------------------------------------------------
% Copyright (c) 2015 Dominic Siedhoff
%-------------------------------------------------------------------------
%
% License: This software may be freely used, shared and modified. It must
%          not be sold. It is provided without any explicit or implicit
%          warranty of any kind. This license text must be included with
%          every copy made.
%
%-------------------------------------------------------------------------




if nargin < 3 || isempty(in_beta)
    in_beta = 1;
end

if nargin < 4 || isempty(in_kDensity)
    in_kDensity = 5;
end

if nargin < 5 || isempty(in_kSMOTE)
    in_kSMOTE = 5;
end

if nargin < 6 || isempty(in_featuresAreNormalized)
    in_featuresAreNormalized = true;
end


if in_beta == 0
    %nothing needs to be done because beta==0 is defined to mean that
    %current class ratio is kept:
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
end


if ~all(in_labels==0 | in_labels==1)
    error('ADASYN: in_labels may contain only the values 0 and 1.');
end


numZeros = sum(in_labels==0);
numOnes  = sum(in_labels==1);

if numOnes == numZeros
    %nothing needs to be done because if classes are already balanced, then
    %for any in_beta, this is already the desired result.
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
else
    if numZeros > numOnes
        majLabel = logical(0);
        minLabel = logical(1);
    else
        majLabel = logical(1);
        minLabel = logical(0);
    end
end

%rename:
S = in_features;
clear in_features;

%feature sets by class:
Smin = S(in_labels==minLabel,:);
Smaj = S(in_labels==majLabel,:);

%eq (3):
G = (size(Smaj,1) - size(Smin,1)) * in_beta;
G = round(G);

%handle boundary cases:
if size(Smin,1)==0
    warning('ADASYN: there were no examples of the minority class in the data. hence balancing is not possible. Returning empty matrices.');
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
end

if size(Smin,1)==1
    warning('ADASYN: there was only one example of the minority class in the data. Hence returning G copies of that single example for balancing.');
    out_featuresSyn = repmat(Smin, [G 1]);
    out_labelsSyn   = logical(minLabel * ones([G 1]));
    return;
end


if in_featuresAreNormalized
    knnDistance = 'euclidean';
else
    %IMPORTANT: using 'seuclidean' as the distance measure means that
    %standardized Euclidean distance is used, i.e. the standard deviation
    %of the coordinates is automatically divided away. hence, using
    %'seuclidean' instead of 'euclidean' saves the effort of normalizing
    %the feature values by a Z-transformation (0mean,1var).
    %cf. documentation of input parameter in_featuresAreNormalized for more
    %information.
    knnDistance = 'seuclidean';
end

%kNN for density estimation:
idcs = knnsearch_nonflat(S,Smin, 'K',in_kDensity+1, 'Distance',knnDistance);
%note: why in_kDensity+1? because Smin is a subset of S and hence all
%points in Smin have a trivial nearest neighbor in S with distance 0.
%but that neighbor is not interesting because it's the point from Smin
%itself. hence remove it:
idcs = idcs(:,2:end);


%compute the \Gamma values (eq. (4) in reference [2]):
Gamma = zeros([size(Smin,1) 1]);
for cmi=1:size(Smin,1)  %cmi: current minority example index
    cNNs = idcs(cmi,:);             %current NearestNeighbors
    cNNsLabels = in_labels(cNNs);   %labels of cNNs:
    cNNsLabelsMaj = (cNNsLabels == majLabel);
    cDelta = nnz(cNNsLabelsMaj);    %the Delta_i of eq. (4) in reference [2]
    
    %write Gamma, not yet normalized:
    Gamma(cmi) = cDelta / in_kDensity;

end

%normalize Gamma to give a distribution function:
if sum(Gamma)==0
    %if there is no class overlap w.r.t. these knn settings, create a
    %uniform distribution:
    Gamma = 1/length(Gamma) * ones(size(Gamma));
else
    %create nonuniform distribution:
    Gamma = Gamma / sum(Gamma);     %this makes it exactly eq. (4) in reference [2]
end

%compute g_i (eq. (5) in reference [2]):
%these g_i are the numbers of synthetic examples to be generated from each
%example in Smin
g = round(Gamma * G);

if sum(g)==0
    warning('ADASYN: Classes are already well-balanced (i.e. sum(g)==0). Returning empty matrices.');
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
end

%with this g known, call the ADASYN_SMOTE subroutine...:
out_featuresSyn = ADASYN_SMOTE(Smin,g,in_kSMOTE,knnDistance);
%...and generate the labels:
out_labelsSyn = logical(minLabel * ones([size(out_featuresSyn,1) 1]));



function Ssyn = ADASYN_SMOTE(Smin,g,k,knnDistance)
%subroutine implementing SMOTE algorithm as it is to be used by function
%ADASYN(). cf. section 3.1.3 in the following paper for details:
%
%[2]: H. He and E.A. Garcia, "Learning from imbalanced data",
%Knowledge and Data Engineering, IEEE Transactions on 21, no. 9,
%pp. 1263--1284, (2009).
%
%
% INPUTS:
%----------
%
%Smin:
%minority set from ADASYN()
%
%g:
%minority example synthesis counts as computed by ADASYN()
%
%k:
%number of neighbors to be regarded in kNN for SMOTE algorithm
%
%knnDistance:
%distance function used in kNN for SMOTE algorithm. depends on ADASYN's
%parameter in_featuresAreNormalized. please type "help ADASYN" into
%MATLAB's console for more information
%
%
% OUTPUTS:
%----------
%Ssyn: set of synthetic examples created from input set Smin by applying
%the SMOTE algorithm


%determine nearest neighbors:
idcs = knnsearch_nonflat(Smin,Smin, 'K',k+1, 'Distance',knnDistance);
%note: why k+1? because we search kNNs of Smin in Smin itself and hence all
%points in Smin have a trivial nearest neighbor in Smin with distance 0.
%but that neighbor is not interesting because it's the point from Smin
%itself. hence remove it:
idcs = idcs(:,2:end);

%initialize output and writing target as an empty matrix
Ssyn = zeros([0 size(Smin,2)]);

%for every minority example xi...
for cei=1:size(Smin,1)  %cei: current example index
    
    %current minority example:
    xi = Smin(cei,:);
    
    %number of synthetic examples to be created from xi:
    gi = g(cei);
    
    %allocate space for gi examples to be created from xi:
    xiSyn = zeros(gi, size(Smin,2));
    
    %...iterate over synthetic examples to be created from xi and
    %random partner from set of nearest neighbors:
    for csi=1:gi    %csi: current synthetic example index
        
        %get random partner example from nearest neighbors of xi:
        %neighbor index:
        nIdx = idcs(cei, randi(size(idcs,2)));
        %neighbor:
        xiHat = Smin(nIdx,:);
        
        %create synthetic example as according to eq. (1) in reference [2]:
        delta = rand(1);
        xSyn = xi + delta * (xiHat - xi);
        
        %write it to xiSyn:
        xiSyn(csi,:) = xSyn;
        
    end
    
    %append examples synthesized from xi to overall synthetic example set:
    Ssyn = [Ssyn; xiSyn];
    
end





function [IDX,D] = knnsearch_nonflat(X,Y, varargin)
%wraps knnsearch from MATLAB's statistics toolbox.
%knnsearch_nonflat executes knnsearch only on the dimensions with nonzero
%standard deviation, i.e. the flat dimensions are not passed on to
%knnsearch. this prevents pdist2 from producing the following warning in
%the context of standardized Euclidean distance ('seuclidean') in the
%presence of flat dimensions:
%"Warning: Some columns of S are zeros."
%if this warning occurs, pdist2 (and as a consequence knnsearch) gives only
%bad dummy results because the standardized Euclidean distance can not be
%computed properly in the presence of flat dimensions.
%using knnsearch_nonflat prevents this by filtering out flat dimensions.

nonflatX = std(X) ~= 0;
nonflatY = std(Y) ~= 0;

nonflat = nonflatX & nonflatY;

[IDX,D] = knnsearch(X(:,nonflat), Y(:,nonflat), varargin{:});


