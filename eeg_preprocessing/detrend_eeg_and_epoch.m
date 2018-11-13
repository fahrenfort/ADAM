function detrend_eeg_and_epoch(filepath,filenames,outpath, polynomial_order, epoch_pad, start_epoch, end_epoch, varargin)
% function detrend_eeg_and_epoch(filepath,filename,outpath, polynomial_order, epoch_pad, start_epoch, end_epoch, varargin)
% reads in EEGLAB data, detrends, epochs and writes out. 
% filenames: names of files either in a cell array or as comma separated
% string. Wildcards * and ? can be used, e.g. filenames = '*' will
% take all the .set files in the input filepath as sources. 
%
%       polynomial order        order of the polynomial that is used for detrending (Default: 30).
%       epoch_pad               used to pad each side of the original epoch, e.g. epoch_pad of 100
%                               creates 200 sec epochs (Default: 100). This padding window is only
%                               used during detrending and is discarded when doing the final
%                               epoching.
%       start_epoch             in seconds
%       end_epoch               in seconds
% 
% varargin is list of conditions to epoch, either as strings or numbers
% separated by comma's, or as a cell array
% Use  cond_string([1001:1024 1101:1124]) when called from create_qsub_files in
% order to epoch all conditions running from 1001 to 1024 and from 1101 to
% 1124. 
%
% example: detrend_eeg_and_epoch('c:\alldata\','subject01', 'c:\filtereddata', 30, 100, -.25, 0, cond_string([1001:1024 1101:1124]));
%
% J.J.Fahrenfort, VU 2018

% some input checking
if nargin < 5 || isempty(polynomial_order)
    polynomial_order = 30;
end
if nargin < 4 || isempty(epoch_pad)
    epoch_pad = 100;
end
if nargin < 3
    error('not enough arguments');
end
if isempty(outpath)
    outpath = filepath;
else
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
end
if ischar(hicutoff);
    hicutoff = string2double(hicutoff);
end
if ischar(locutoff)
    locutoff = string2double(locutoff);
end
if ischar(start_epoch);
    start_epoch = string2double(start_epoch);
end
if ischar(end_epoch)
    end_epoch = string2double(end_epoch);
end
if ~iscell(filenames) && (~isempty(strfind(filenames,'*')) || ~isempty(strfind(filenames,'?')))
    if ~strcmp(filenames(end-3:end),'.set')
        filenames = [filenames '.set'];
    end
    filenames = dir([filepath filesep filenames]);
    filenames = {filenames(:).name};
end
if ~iscell(filenames)
    filenames = regexp(filenames, ',', 'split');
end
% conditions can either be defined as a string of comma separated values,
% or a cell array of numbers or strings, and is converted to a ell array of
% strings to be fed into pop_epoch
if numel(varargin) == 1
    if ~iscell(varargin{1})
        conditions = regexp(varargin{1}, ',', 'split');
    end
else
    conditions = varargin;
    if ~isempty(conditions)
        if ~ischar(conditions{1})
            conditions = regexp(vec2str([conditions{:}]), ',', 'split');
        end
    end
end
% go
for filename = filenames
    [~,fname,~] = fileparts(filename{1});
    % load and filter
    EEG = pop_loadset('filename',[fname '.set'],'filepath',filepath);
    
    % mirror-pad edges of the unepoched data, so that extracting full trials will not be problematic
    EEG.data = padarray(EEG.data,[0 epoch_pad*EEG.srate],'both','symmetric');
    tdiff = EEG.times(2)-EEG.times(1);
    EEG.times(end+1:end+2*epoch_pad*EEG.srate) = EEG.times(end)+tdiff:tdiff:EEG.times(end)+tdiff*(2*epoch_pad*EEG.srate);
    EEG.pnts = length(EEG.times);
    EEG.xmin = EEG.times(1);
    EEG.xmax = EEG.times(end);
    % also fix latencies after padding
    for evi=1:length(EEG.event)
        EEG.event(evi).latency = EEG.event(evi).latency+EEG.srate*epoch_pad;
        EEG.urevent(evi).latency = EEG.urevent(evi).latency+EEG.srate*epoch_pad;
    end

    % Wide epoch to perform the detrending on
    detrendtime     = [start_epoch end_epoch] + [-1*epoch_pad epoch_pad];
    EEG = pop_epoch( EEG, conditions, detrendtime);
    
    
    tic
    fprintf('Robust detrending @ %ith order polynomial over an interval of [%s] sec.\n',polynomial_order,vec2str(detrendtime));
    
    all_events = pop_exportepoch(EEG);
    
    filt_dat = zeros(size(EEG.data)); % EEG.data contains channels x time x epochs
    nChannels = size(EEG.data,1);
    nTimes = size(EEG.data,2);
    nTrials = size(EEG.data,3);
    for cTrial = 1:nTrials
        
        wt = ones(nTimes,1); % create a mask
        
        % remove the entire trial from the mask, running from 0 until 250 ms prior to the next trial
        % (effectively catching the next baseline period)
        if cTrial == nTrials
            timeOfNextEvent = end_epoch*1000 + 500; % for the last trial just go to the end epoch (multiply by 1000 to obtain ms) and add 500 ms
        else
            indNextEvent = find(EEG.epoch(cTrial).event,all_events(cTrial+1));
            timeOfNextEvent = EEG.epoch(cTrial).eventlatency(indNextEvent) - 250;
        end
        tInd = nearest(EEG.times,[0 timeOfNextEvent]); 
        wt(tInd(1):tInd(2))=0; % mask that interval out by setting it to zero 
        
        %[tmp,w1] = nt_detrend(EEG.data(:,:,cTrial)',1,repmat(wt,[1 nTrials])); % start with 1st order
        % Why start with first order?? Just go right away!
        if polynom>1
            [tmp, ~] = nt_detrend(tmp,polynomial_order,w1); % then nth order with mask of previous step
            filt_dat(:,:,cTrial) = tmp';
        else
            filt_dat(:,:,cTrial) = tmp';
        end
    end
    toc

    EEG.data = filt_dat;
    clear filt_dat
    
    % narrow epoch
    EEG = pop_select( EEG, 'time',[start_epoch end_epoch]);
    
    % OLD pop_epoch:
    % EEG = pop_epoch(EEG, conditions, [start_epoch end_epoch], 'newname', ['epoched_ ' fname], 'epochinfo', 'yes');

    % save data
    pop_saveset(EEG, 'filename',[fname '.set'],'filepath',outpath);
end


