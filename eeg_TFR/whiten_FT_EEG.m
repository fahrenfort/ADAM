function [FT_EEG, IE] = whiten_FT_EEG(FT_EEG,condSet,IE,field)
% whiten_FT_EEG whitens the EEG data contained in an FT_EEG struct using Mahalanobis or ZCA
% whitening, see https://en.wikipedia.org/wiki/Whitening_transformation. This is also called
% sphering, or Multivariate Noise Normalization (MNN) also see:
% https://www.biorxiv.org/content/early/2017/08/06/172619
%
% First computes covariance matrix per time point and per stimulus class. Averages covariance
% matrices using across time points and across stimulus classes and inverts this matrix. Next, it
% uses the inverted matrix to whiten the original data, so that the covariance is the identity
% matrix, a set of uncorrelated variables with unit variance. An inverted matrix can also be
% supplied when calling this function, in this case the data will be whitened using the supplied
% matrix.
%
% Part of the ADAM toolbox, by J.J. Fahrenfort, 2018
%
% see also: covCor

if nargin < 4
    field = 'trial';
end
if nargin < 3
    IE = [];
end

if ~isfield(FT_EEG,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on FT_EEG to resolve this.');
else
    dims = regexp(FT_EEG.dimord, '_', 'split');
    chandim = find(strcmp(dims,'chan'));
    timedim = find(strcmp(dims,'time'));
    trialdim = find(strcmp(dims,'rpt'));
end

if isempty(chandim | timedim | trialdim) || numel(dims) > 3
    error('incorrect dimensions: should have time, channel and trial');
end
if ~isfield(FT_EEG,field)
    error([field ' should contain the data, but ' field ' is not a field of FT_EEG']);
end

% get trial info
trialinfo = FT_EEG.trialinfo;
% data
data = permute(FT_EEG.(field),[trialdim chandim timedim]);
[~, nChannels, nTime] = size(data);
nConds = numel(condSet);

% compute IE (if it was not passed to the function)
if isempty(IE)
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
    for cTrial = 1:numel(trialinfo)
        X=squeeze(data(cTrial,:,cT)); % for each trial and time point
        W=X*IE; % WHITEN
        data_new(cTrial,:,cT)=W; % whitened across channels, for each condition and time point
    end
end

FT_EEG.dimord = 'rpt_chan_time'; % should be 'rpt_chan_time' for consistency with conversion to raw format in ft_checkdata
FT_EEG.(field) = data_new; % whitened data