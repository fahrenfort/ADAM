function [FT_EEG, IE] = whiten_FT_EEG(FT_EEG,condSet,IE,field,avgtime)
% whiten_FT_EEG whitens the EEG data contained in an FT_EEG struct using Mahalanobis or ZCA
% whitening, see https://en.wikipedia.org/wiki/Whitening_transformation. This is also called
% sphering, or Multivariate Noise Normalization (MNN) also see:
% https://www.biorxiv.org/content/early/2017/08/06/172619
%
% First computes covariance matrix per time point and per stimulus class. Averages covariance
% matrices across stimulus classes and optionally across time points, and inverts this matrix. Next,
% it uses the inverted matrix to whiten the original data, so that the covariance is the identity
% matrix, a set of uncorrelated variables with unit variance. An inverted matrix can also be
% supplied when calling this function, in this case the data will be whitened using the supplied
% matrix.
%
% Internal function of the ADAM toolbox, by J.J.Fahrenfort, 2018
%
% See also: ADAM_MVPA_FIRSTLEVEL, COVCOR

if nargin < 5
    avgtime = false;
end
if nargin < 4
    field = 'trial';
end
if nargin < 3
    IE = [];
end

if ~isfield(FT_EEG,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on FT_EEG to resolve this.');
end
if ~isfield(FT_EEG,field)
    error([field ' should contain the data, but ' field ' is not a field of FT_EEG']);
end
if isempty(IE)
    computeIE = true;
else
    computeIE = false;
end

% get trial info
trialinfo = FT_EEG.trialinfo;

% dimord should be 'rpt_chan_time' in this function, could rewrite to make sure input format does not change to save memory
FT_EEG = fix_dimord(FT_EEG,'rpt_chan_time');
data = FT_EEG.(field);
[~, nChannels, nTime] = size(data);
nConds = numel(condSet);

% compute IE (if it was not passed to the function)
if avgtime && computeIE
    % pre-allocate (time, conditions, channels, channels)
    E=NaN(nTime,nConds,nChannels,nChannels);
    for cT = 1:nTime
        for cCondSet = 1:nConds
            thisCondSet = [condSet{cCondSet}];
            if ischar(thisCondSet)
                thisCondSet = string2double(thisCondSet);
            end
            X = squeeze(data(ismember(trialinfo,thisCondSet),:,cT)); % X = trials x channels
            E(cT,cCondSet,:,:)=covCor(X); % E = time x condition/class x channels x channels
        end
    end
    EE = squeeze(mean(mean(E,1),2)); % average over time and condition/class
    clear E;
    % invert E -> this is Mahalanobis or ZCA whitening, see https://en.wikipedia.org/wiki/Whitening_transformation
    IE=EE^(-0.5);  % IE is channels x channels
end
% figure; imagesc(IE); % plot if you are interested

% pre-allocate memory (trials, channels, time)
data_new=NaN(numel(trialinfo),nChannels,nTime);
for cT = 1:nTime
    if ~avgtime && isempty(IE)
        E=NaN(nConds,nChannels,nChannels);
        for cCondSet = 1:nConds
            thisCondSet = [condSet{cCondSet}];
            if ischar(thisCondSet)
                thisCondSet = string2double(thisCondSet);
            end
            X = squeeze(data(ismember(trialinfo,thisCondSet),:,cT)); % X = trials x channels
            E(cCondSet,:,:)=covCor(X); % E = condition/class x channels x channels
        end
        EE = squeeze(mean(E,1)); % only average over condition/class
        IE=EE^(-0.5);  % IE is channels x channels
    end
    for cTrial = 1:numel(trialinfo)
        X=squeeze(data(cTrial,:,cT)); % for each trial and time point
        W=X*IE; % WHITEN
        data_new(cTrial,:,cT)=W; % whitened across channels, for each condition and time point
    end
end
FT_EEG.(field) = data_new; % whitened data