function subtract_conditions_TFR(filepath,filename,outpath,varargin)
% function FT_EEG = subtract_conditions_TFR(filepath,filename,outpath,condSet)
% SubtractS conditions in condSet inside FT_EEG. 
% filename contains a single filename
% outpath contains the folder where results should be stored (if empty
% defaults to filepath)condSet can be set up as follows (example): 
% varargin is a variable set of conditions to select from, by which the
% main conditions that should be subtracted are specified as string with
% comma separated values: 
% cond1 = '1,2';
% cond2 = '3,4';
% Subtracts condition 2 from condition 1, becoming the 'new' condition 1
% and condition 4 from condition 3, becoming the 'new' condition 2.
% example usage:
% subtract_conditions_TFR('/Users/VU-MBP/Documents/TFR_TEMP','TFR_BIN_TOT_subj02_AB.mat','','1,2','3,4','5,6','7,8','9,10');
%
% J.J.Fahrenfort, VU 2015

if nargin < 4
    error('this function requires at least 4 arguments')
end

% some bookkeeping
filename = regexp(filename, ',', 'split');
if numel(filename) > 1
    error('multiple files are not yet implemented');
else
    filename = filename{1};
    filename(strfind(filename,'.mat'):end) = [];
end
if isempty(outpath)
    outpath = filepath;
else
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
end
if ~iscell(varargin{1})
    for cCond = 1:numel(varargin)
        condSet{cCond} = str2double(varargin{cCond});
    end
else
    condSet = varargin{:};
end

% load and check data
TFR = load([filepath filesep filename '.mat']);
if ~isfield(TFR,'dimord')
    error('the input dataset is not in the required standard fieldtrip format. you might want to run ft_timelockbaseline on TFR to resolve this.');
else
    dims = regexp(TFR.dimord, '_', 'split');
    chandim = find(strcmp(dims,'chan'));
    timedim = find(strcmp(dims,'time'));
    trialdim = find(strcmp(dims,'rpt'));
    freqdim = find(strcmp(dims,'freq'));
end

if isempty(chandim | timedim | trialdim | isempty(freqdim)) || numel(dims) > 4
    error('incorrect dimensions: should have time, channel and trial, cannot subtract TFR');
end

% get the data
oldpowspctrm = permute(TFR.powspctrm,[chandim timedim trialdim freqdim]);
newdimord = 'chan_time_rpt_freq';
oldtrialinfo =  TFR.trialinfo;

% find indexes of conditions and subtract

newpowspctrm = [];%zeros(size(oldpowspctrm));
newtrialinfo = [];
for cCond = 1:numel(condSet)
    cond2Subtract = condSet{cCond};
    if numel(cond2Subtract) ~= 2
        error('can only subtract if exactly two conditions are present');
    end
    cond1Index = find(oldtrialinfo == cond2Subtract(1));
    cond2Index = find(oldtrialinfo == cond2Subtract(2));
    maxSize = min([numel(cond1Index) numel(cond2Index)]);
    cond1Index = cond1Index(1:maxSize);
    cond2Index = cond2Index(1:maxSize);
    newpowspctrm(:,:,end+1:end+maxSize,:) = oldpowspctrm(:,:,cond1Index,:)-oldpowspctrm(:,:,cond2Index,:);
    newtrialinfo(end+1:end+maxSize,1) = cCond;
end

% replace relevant field
TFR.dimord = newdimord;
TFR.powspctrm = newpowspctrm;
TFR.trialinfo = newtrialinfo;
if isfield(TFR,'oldindex')
    TFR = rmfield(TFR,'oldindex');
end

% save outcome
save([outpath filesep 'SUBTR_' filename '.mat'],'-v7.3','-struct','TFR');

