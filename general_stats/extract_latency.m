function latencies = extract_latency(cfg, stats)
% internal function of the ADAM toolbox that uses the function latency.m to
% compute GA, Jackmiller and individual latencies.

% first get the values to operate on
accs = stats.indivClassOverTime;                     % classifier performance or muV data
% now subtract the chance level, or the latency function will not work properly
if isfield(cfg,'chance')
    accs = accs - cfg.chance;
end
% we should fix this if not operating on the diagonal
if ~ismatrix(accs)
    error('Need to introduce some code to handle temporal generalization matrices and/or time-frequency matrices when computing onset.');
end
accs = reshape(accs,[size(accs,1) 1 size(accs,2)]);     % subjects * electrode * time
times = stats.settings.times{1};                        % in seconds

% compute some defaults
if ~isfield(cfg,'peakWin')
    % find first peak
    if ~isempty(stats.pStruct.posclusters) && ~isempty(stats.pStruct.negclusters)
        if stats.pStruct.negclusters(1).peak_time < stats.pStruct(1).posclusters.peak_time
            peakfield = 'negclusters';
        else
            peakfield = 'posclusters';
        end
    elseif ~isempty(stats.pStruct.posclusters)
        peakfield = 'posclusters';
    elseif ~isempty(stats.pStruct.negclusters)
        peakfield = 'negclusters';
    else
        latencies = [];
        return
    end
    if strcmpi(peakfield,'negclusters')
        sign = -1;
        peakWin = [stats.pStruct.negclusters(1).start_time stats.pStruct.negclusters(1).stop_time ]; % signficant window
    elseif strcmpi(peakfield,'posclusters')
        sign = 1;
        peakWin = [stats.pStruct.posclusters(1).start_time stats.pStruct.posclusters(1).stop_time ]; % signficant window
    end
end
chans = 1;                                          % only one channel
extract = 'onset';                                  % percentage amplitude reached
percAmp = .5;                                       % 50% amplitude
fig = false;                                        % set to true to plot figures
%sampRate = (numel(times)-1)/(times(end)-times(1));  % in Hz
times = times * 1000;                               % now in milliseconds
peakWidth = 10;                                     % average across 10 ms. to find peak
warnings = false;

% unpack settings
v2struct(cfg);

% insert true values (after defaults have been set)
cfg = [];                                               % reset config
cfg.chans = chans;                                      % only one channel
cfg.extract = extract;                                  % percentage amplitude reached
cfg.percAmp = percAmp;                                  % 50% amplitude
cfg.sign = sign;                                        % is peak positive or negative?
cfg.fig = fig;                                          % set to true to plot figures
cfg.times = times;                                      % in milliseconds
cfg.timeFormat = 'ms';
%cfg.sampRate = sampRate;                                % in Hz
cfg.peakWin = peakWin;                                  % cluster based permutation window
cfg.peakWidth = peakWidth;                              % in milliseconds
cfg.warnings = warnings;

% compute relevant values
cfg.aggregate = 'GA';                                   % do Grand Average
latencies.GA = latency(cfg,accs);
cfg.aggregate = 'individual';                           % do individual onsets
latencies.individual = latency(cfg,accs);
cfg.aggregate = 'jackMiller';                           % do a jackknife
latencies.jackknife = latency(cfg,accs);

