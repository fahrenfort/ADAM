function [cfg, artifact] = ft_artifact_zvalue_jjf(cfg, data)

% FT_ARTIFACT_ZVALUE reads the interesting segments of data from file and
% identifies artifacts by means of thresholding the z-transformed value
% of the preprocessed raw data. Depending on the preprocessing options,
% this method will be sensitive to EOG, muscle or jump artifacts.
% This procedure only works on continuously recorded data.
%
% Use as
%   [cfg, artifact] = ft_artifact_zvalue(cfg)
% or
%   [cfg, artifact] = ft_artifact_zvalue(cfg, data)
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% If you are calling FT_ARTIFACT_ZVALUE with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%   cfg.headerfile
%   cfg.datafile
% and optionally
%   cfg.headerformat
%   cfg.dataformat
%
% If you are calling FT_ARTIFACT_ZVALUE with also the second input argument
% "data", then that should contain data that was already read from file
% a call to FT_PREPROCESSING.
%
% If you encounter difficulties with memory usage, you can use
%   cfg.memory = 'low' or 'high', whether to be memory or computationally efficient, respectively (default = 'high')
%
% The required configuration settings are:
%   cfg.trl
%   cfg.continuous
%   cfg.artfctdef.zvalue.channel
%   cfg.artfctdef.zvalue.cutoff
%   cfg.artfctdef.zvalue.trlpadding
%   cfg.artfctdef.zvalue.fltpadding
%   cfg.artfctdef.zvalue.artpadding
%
% If you specify
%   cfg.artfctdef.zvalue.interactive = 'yes', a GUI will be started and you
%     can manually accept/reject detected artifacts, and/or change the threshold
%
% If you specify
%   cfg.artfctdef.zvalue.artfctpeak='yes', the maximum value of the artifact
%       within its range will be found; saved into cfg.artfctdef.zvalue.peaks
% Use also, e.g. as input to DSS option of ft_componentanalysis:
%   cfg.artfctdef.zvalue.artfctpeakrange= [-0.25 0.25]; (in seconds), (for example)
%       to indicate range around peak to include, saved into cfg.artfctdef.zvalue.dssartifact
%       Default values: [0 0]
%       Range will respect trial boundaries (i.e. be shorter if peak is near beginning or end of trial)
%       Samples between trials will be removed; thus this won't match .sampleinfo of the data structure
%
% Configuration settings related to the preprocessing of the data are
%   cfg.artfctdef.zvalue.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.artfctdef.zvalue.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.artfctdef.zvalue.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.artfctdef.zvalue.lnfilter      = 'no' or 'yes'  line noise removal using notch filter
%   cfg.artfctdef.zvalue.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.artfctdef.zvalue.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.artfctdef.zvalue.lpfreq        = lowpass  frequency in Hz
%   cfg.artfctdef.zvalue.hpfreq        = highpass frequency in Hz
%   cfg.artfctdef.zvalue.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.artfctdef.zvalue.lnfreq        = line noise frequency in Hz, default 50Hz
%   cfg.artfctdef.zvalue.lpfiltord     = lowpass  filter order
%   cfg.artfctdef.zvalue.hpfiltord     = highpass filter order
%   cfg.artfctdef.zvalue.bpfiltord     = bandpass filter order
%   cfg.artfctdef.zvalue.lnfiltord     = line noise notch filter order
%   cfg.artfctdef.zvalue.medianfiltord = length of median filter
%   cfg.artfctdef.zvalue.lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.artfctdef.zvalue.hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.artfctdef.zvalue.bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.artfctdef.zvalue.detrend       = 'no' or 'yes'
%   cfg.artfctdef.zvalue.demean        = 'no' or 'yes'
%   cfg.artfctdef.zvalue.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.artfctdef.zvalue.hilbert       = 'no' or 'yes'
%   cfg.artfctdef.zvalue.rectify       = 'no' or 'yes'
%
% See also FT_REJECTARTIFACT, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG, FT_ARTIFACT_EOG,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_ZVALUE

% Copyright (C) 2003-2011, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_artifact_zvalue.m 10259 2015-02-24 16:03:12Z roboos $

revision = '$Id: ft_artifact_zvalue.m 10259 2015-02-24 16:03:12Z roboos $';

% % do the general setup of the function
% ft_defaults
% ft_preamble init
% ft_preamble provenance
% ft_preamble loadvar data
% 
% % the abort variable is set to true or false in ft_preamble_init
% if abort
%   return
% end

% set default rejection parameters
cfg.headerformat                 = ft_getopt(cfg,                  'headerformat', []);
cfg.dataformat                   = ft_getopt(cfg,                  'dataformat',   []);
cfg.memory                       = ft_getopt(cfg,                  'memory',       'high');
cfg.artfctdef                    = ft_getopt(cfg,                  'artfctdef',    []);
cfg.artfctdef.zvalue             = ft_getopt(cfg.artfctdef,        'zvalue',       []);
cfg.artfctdef.zvalue.method      = ft_getopt(cfg.artfctdef.zvalue, 'method',       'all');
cfg.artfctdef.zvalue.ntrial      = ft_getopt(cfg.artfctdef.zvalue, 'ntrial',       10);
cfg.artfctdef.zvalue.channel     = ft_getopt(cfg.artfctdef.zvalue, 'channel',      {});
cfg.artfctdef.zvalue.trlpadding  = ft_getopt(cfg.artfctdef.zvalue, 'trlpadding',   0);
cfg.artfctdef.zvalue.fltpadding  = ft_getopt(cfg.artfctdef.zvalue, 'fltpadding',   0);
cfg.artfctdef.zvalue.artpadding  = ft_getopt(cfg.artfctdef.zvalue, 'artpadding',   0);
cfg.artfctdef.zvalue.interactive = ft_getopt(cfg.artfctdef.zvalue, 'interactive',  'no');
cfg.artfctdef.zvalue.cumulative  = ft_getopt(cfg.artfctdef.zvalue, 'cumulative',   'yes');
cfg.artfctdef.zvalue.artfctpeak  = ft_getopt(cfg.artfctdef.zvalue, 'artfctpeak',   'no');
cfg.artfctdef.zvalue.artfctpeakrange  = ft_getopt(cfg.artfctdef.zvalue, 'artfctpeakrange',[0 0]);

% for backward compatibility
cfg.artfctdef        = ft_checkconfig(cfg.artfctdef,        'renamed', {'blc',      'demean'});
cfg.artfctdef        = ft_checkconfig(cfg.artfctdef,        'renamed', {'blcwindow' 'baselinewindow'});
cfg.artfctdef.zvalue = ft_checkconfig(cfg.artfctdef.zvalue, 'renamed', {'sgn',      'channel'});
cfg.artfctdef.zvalue = ft_checkconfig(cfg.artfctdef.zvalue, 'renamed', {'feedback', 'interactive'});

if isfield(cfg.artfctdef.zvalue, 'artifact')
  fprintf('zvalue artifact detection has already been done, retaining artifacts\n');
  artifact = cfg.artfctdef.zvalue.artifact;
  return
end

% set feedback
cfg.feedback = ft_getopt(cfg, 'feedback',   'text');

% flag whether to compute z-value per trial or not, rationale being that if
% there are fluctuations in the variance across trials (e.g. due to
% position differences in MEG measurements) which don't have to do with the artifact per se,
% the detection is compromised (although the data quality is questionable
% when there is a lot of movement to begin with).
pertrial    = strcmp(cfg.artfctdef.zvalue.method, 'trial');
demeantrial = strcmp(cfg.artfctdef.zvalue.method, 'trialdemean');
if pertrial
  if isfield(cfg.artfctdef.zvalue, 'ntrial') && cfg.artfctdef.zvalue.ntrial>0
    pertrial = cfg.artfctdef.zvalue.ntrial;
  else
    error('you should specify cfg.artfctdef.zvalue.ntrial, and it should be > 0');
  end
end

% the data is either passed into the function by the user or read from file with cfg.inputfile
hasdata = exist('data', 'var');

if ~hasdata
  % only cfg given, read data from disk
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
  hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat);
  trl = cfg.trl;
  
else
  % check whether the value for trlpadding makes sense
  % negative trlpadding only allowed with in-memory data
  %if cfg.artfctdef.zvalue.trlpadding < 0
  %  error('negative trlpadding is only allowed with in-memory data');
  %end
  % check if the input data is valid for this function
  data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');
  cfg = ft_checkconfig(cfg, 'forbidden', {'dataset', 'headerfile', 'datafile'});
  hdr = ft_fetch_header(data);
  trl = data.sampleinfo;
end

% set default cfg.continuous
if ~isfield(cfg, 'continuous')
  if hdr.nTrials==1
    cfg.continuous = 'yes';
  else
    cfg.continuous = 'no';
  end
end

trlpadding    = round(cfg.artfctdef.zvalue.trlpadding*hdr.Fs);
fltpadding    = round(cfg.artfctdef.zvalue.fltpadding*hdr.Fs);
artpadding    = round(cfg.artfctdef.zvalue.artpadding*hdr.Fs);
trl(:,1)      = trl(:,1) - trlpadding;       % pad the trial with some samples, in order to detect
trl(:,2)      = trl(:,2) + trlpadding;       % artifacts at the edges of the relevant trials.
if size(trl, 2) >= 3
  trl(:,3)      = trl(:,3) - trlpadding;     % the offset can ofcourse be adjusted as well
elseif hasdata
  % reconstruct offset
  for tr=1:size(trl, 1)
    % account for 0 might not be in data.time
    t0         = interp1(data.time{tr}, 1:numel(data.time{tr}), 0, 'linear', 'extrap');
    trl(tr, 3) = -t0+1 - trlpadding;
  end
else
  % assuming that the trial starts at t=0s
  trl(:, 3) = trl(:, 1);
end
trllength     = trl(:,2) - trl(:,1) + 1;     % length of each trial
numtrl        = size(trl,1);
cfg.artfctdef.zvalue.trl = trl;              % remember where we are going to look for artifacts
cfg.artfctdef.zvalue.channel = ft_channelselection(cfg.artfctdef.zvalue.channel, hdr.label);
sgnind        = match_str(hdr.label, cfg.artfctdef.zvalue.channel);
numsgn        = length(sgnind);
thresholdsum  = strcmp(cfg.artfctdef.zvalue.cumulative, 'yes');

if numsgn<1
  error('no channels selected');
end

% read the data and apply preprocessing options
sumval = zeros(numsgn, 1);
sumsqr = zeros(numsgn, 1);
numsmp = zeros(numsgn, 1);
ft_progress('init', cfg.feedback, ['searching for artifacts in ' num2str(numsgn) ' channels']);
for trlop = 1:numtrl
  ft_progress(trlop/numtrl, 'searching in trial %d from %d\n', trlop, numtrl);
  
  if strcmp(cfg.memory, 'low') % store nothing in memory
    if hasdata
      dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'), 'skipcheckdata', 1);
    else
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'), 'dataformat', cfg.dataformat);
    end
    dat = preproc(dat, cfg.artfctdef.zvalue.channel, offset2time(0, hdr.Fs, size(dat,2)), cfg.artfctdef.zvalue, fltpadding, fltpadding);
    
    if trlop==1 && ~pertrial
      sumval = zeros(size(dat,1), 1);
      sumsqr = zeros(size(dat,1), 1);
      numsmp = zeros(size(dat,1), 1);
      numsgn = size(dat,1);
    elseif trlop==1 && pertrial
      sumval = zeros(size(dat,1), numtrl);
      sumsqr = zeros(size(dat,1), numtrl);
      numsmp = zeros(size(dat,1), numtrl);
      numsgn = size(dat,1);
    end
    
    if ~pertrial
      % accumulate the sum and the sum-of-squares
      sumval = sumval + sum(dat,2);
      sumsqr = sumsqr + sum(dat.^2,2);
      numsmp = numsmp + size(dat,2);
    else
      % store per trial the sum and the sum-of-squares
      sumval(:,trlop) = sum(dat,2);
      sumsqr(:,trlop) = sum(dat.^2,2);
      numsmp(:,trlop) = size(dat,2);
    end
  else % store all data in memory, saves computation time
    if hasdata
      dat{trlop} = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'), 'skipcheckdata', 1);
    else
      dat{trlop} = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'), 'dataformat', cfg.dataformat);
    end
    dat{trlop} = preproc(dat{trlop}, cfg.artfctdef.zvalue.channel, offset2time(0, hdr.Fs, size(dat{trlop},2)), cfg.artfctdef.zvalue, fltpadding, fltpadding);
    
    if trlop==1 && ~pertrial
      sumval = zeros(size(dat{1},1), 1);
      sumsqr = zeros(size(dat{1},1), 1);
      numsmp = zeros(size(dat{1},1), 1);
      numsgn = size(dat{1},1);
    elseif trlop==1 && pertrial
      sumval = zeros(size(dat{1},1), numtrl);
      sumsqr = zeros(size(dat{1},1), numtrl);
      numsmp = zeros(size(dat{1},1), numtrl);
      numsgn = size(dat{1},1);
    end
    
    if ~pertrial
      % accumulate the sum and the sum-of-squares
      sumval = sumval + sum(dat{trlop},2);
      sumsqr = sumsqr + sum(dat{trlop}.^2,2);
      numsmp = numsmp + size(dat{trlop},2);
    else
      % store per trial the sum and the sum-of-squares
      sumval(:,trlop) = sum(dat{trlop},2);
      sumsqr(:,trlop) = sum(dat{trlop}.^2,2);
      numsmp(:,trlop) = size(dat{trlop},2);
    end
  end
end % for trlop
ft_progress('close');

if pertrial>1
  sumval = ft_preproc_smooth(sumval, pertrial)*pertrial;
  sumsqr = ft_preproc_smooth(sumsqr, pertrial)*pertrial;
  numsmp = ft_preproc_smooth(numsmp, pertrial)*pertrial;
end

% compute the average and the standard deviation
datavg = sumval./numsmp;
datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);

if strcmp(cfg.memory, 'low')
  fprintf('\n');
end

zmax = cell(1, numtrl);
zsum = cell(1, numtrl);
zindx = cell(1, numtrl);

% create a vector that indexes the trials, or is all 1, in order
% to a per trial z-scoring, or use a static std and mean (used in lines 317
% and 328)
if pertrial
  indvec = 1:numtrl;
else
  indvec = ones(1,numtrl);
end
for trlop = 1:numtrl
  if strcmp(cfg.memory, 'low') % store nothing in memory (note that we need to preproc AGAIN... *yawn*
    fprintf('.');
    if hasdata
      dat = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'));
    else
      dat = ft_read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'), 'dataformat', cfg.dataformat);
    end
    dat = preproc(dat, cfg.artfctdef.zvalue.channel, offset2time(0, hdr.Fs, size(dat,2)), cfg.artfctdef.zvalue, fltpadding, fltpadding);
    zmax{trlop}  = -inf + zeros(1,size(dat,2));
    zsum{trlop}  = zeros(1,size(dat,2));
    zindx{trlop} = zeros(1,size(dat,2));
    
    nsmp          = size(dat,2);
    zdata         = (dat - datavg(:,indvec(trlop)*ones(1,nsmp)))./datstd(:,indvec(trlop)*ones(1,nsmp));  % convert the filtered data to z-values
    zsum{trlop}   = nansum(zdata,1);                   % accumulate the z-values over channels
    [zmax{trlop},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
    zindx{trlop}      = sgnind(ind);                % also remember the channel number that has the largest z-value
  else
    % initialize some matrices
    zmax{trlop}  = -inf + zeros(1,size(dat{trlop},2));
    zsum{trlop}  = zeros(1,size(dat{trlop},2));
    zindx{trlop} = zeros(1,size(dat{trlop},2));
    
    nsmp          = size(dat{trlop},2);
    zdata         = (dat{trlop} - datavg(:,indvec(trlop)*ones(1,nsmp)))./datstd(:,indvec(trlop)*ones(1,nsmp));  % convert the filtered data to z-values
    zsum{trlop}   = nansum(zdata,1);                   % accumulate the z-values over channels
    [zmax{trlop},ind] = max(zdata,[],1);            % find the maximum z-value and remember it
    zindx{trlop}      = sgnind(ind);                % also remember the channel number that has the largest z-value
  end
  % This alternative code does the same, but it is much slower
  %   for i=1:size(zmax{trlop},2)
  %       if zdata{trlop}(i)>zmax{trlop}(i)
  %         % update the maximum value and channel index
  %         zmax{trlop}(i)  = zdata{trlop}(i);
  %         zindx{trlop}(i) = sgnind(sgnlop);
  %       end
  %     end
end % for trlop

if demeantrial
  for trlop = 1:numtrl
    zmax{trlop} = zmax{trlop}-mean(zmax{trlop},2);
    zsum{trlop} = zsum{trlop}-mean(zsum{trlop},2);
  end
end
%for sgnlop=1:numsgn
%  % read the data and apply preprocessing options
%  sumval = 0;
%  sumsqr = 0;
%  numsmp = 0;
%  fprintf('searching channel %s ', cfg.artfctdef.zvalue.channel{sgnlop});
%  for trlop = 1:numtrl
%    fprintf('.');
%    if hasdata
%      dat{trlop} = ft_fetch_data(data,        'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind(sgnlop), 'checkboundary', strcmp(cfg.continuous,'no'));
%    else
%      dat{trlop} = read_data(cfg.datafile, 'header', hdr, 'begsample', trl(trlop,1)-fltpadding, 'endsample', trl(trlop,2)+fltpadding, 'chanindx', sgnind(sgnlop), 'checkboundary', strcmp(cfg.continuous,'no'));
%    end
%    dat{trlop} = preproc(dat{trlop}, cfg.artfctdef.zvalue.channel(sgnlop), hdr.Fs, cfg.artfctdef.zvalue, [], fltpadding, fltpadding);
%    % accumulate the sum and the sum-of-squares
%    sumval = sumval + sum(dat{trlop},2);
%    sumsqr = sumsqr + sum(dat{trlop}.^2,2);
%    numsmp = numsmp + size(dat{trlop},2);
%  end % for trlop
%
%  % compute the average and the standard deviation
%  datavg = sumval./numsmp;
%  datstd = sqrt(sumsqr./numsmp - (sumval./numsmp).^2);
%
%  for trlop = 1:numtrl
%    if sgnlop==1
%      % initialize some matrices
%      zdata{trlop} = zeros(size(dat{trlop}));
%      zmax{trlop}  = -inf + zeros(size(dat{trlop}));
%      zsum{trlop}  = zeros(size(dat{trlop}));
%      zindx{trlop} = zeros(size(dat{trlop}));
%    end
%    zdata{trlop}  = (dat{trlop} - datavg)./datstd;              % convert the filtered data to z-values
%    zsum{trlop}   = zsum{trlop} + zdata{trlop};                 % accumulate the z-values over channels
%    zmax{trlop}   = max(zmax{trlop}, zdata{trlop});             % find the maximum z-value and remember it
%    zindx{trlop}(zmax{trlop}==zdata{trlop}) = sgnind(sgnlop);   % also remember the channel number that has the largest z-value
%
%    % This alternative code does the same, but it is much slower
%    %   for i=1:size(zmax{trlop},2)
%    %       if zdata{trlop}(i)>zmax{trlop}(i)
%    %         % update the maximum value and channel index
%    %         zmax{trlop}(i)  = zdata{trlop}(i);
%    %         zindx{trlop}(i) = sgnind(sgnlop);
%    %       end
%    %     end
%  end
%  fprintf('\n');
%end % for sgnlop

for trlop = 1:numtrl
  zsum{trlop} = zsum{trlop} ./ sqrt(numsgn);
end

% always create figure
h = figure;
set(h, 'visible', 'off');

opt.artcfg       = cfg.artfctdef.zvalue;
opt.artval       = {};
opt.artpadding   = artpadding;
opt.cfg          = cfg;
opt.channel      = 'artifact';
opt.hdr          = hdr;
opt.numtrl       = size(trl,1);
opt.quit         = 0;
opt.threshold    = cfg.artfctdef.zvalue.cutoff;
% ADDED BY JOHANNES TO THRESHOLD AT abs(min(Z)+cutoff) 
if opt.threshold < 0
    if thresholdsum 
        %opt.threshold = median([zsum{:}]) + abs(min([zsum{:}])-median([zsum{:}])) + abs(cfg.artfctdef.zvalue.cutoff);
        opt.threshold = mean([zsum{:}]) + abs(min([zsum{:}])) + abs(cfg.artfctdef.zvalue.cutoff);
    else 
        %opt.threshold = median([zmax{:}]) + abs(min([zmax{:}])-median([zmax{:}])) + abs(cfg.artfctdef.zvalue.cutoff);
        opt.threshold = mean([zmax{:}]) + abs(min([zmax{:}])-median([zmax{:}])) + abs(cfg.artfctdef.zvalue.cutoff);
    end
end
opt.thresholdsum = thresholdsum;
opt.trialok      = true(1,opt.numtrl); % OK by means of objective criterion
opt.keep         = zeros(1,opt.numtrl); % OK overruled by user +1 to keep, -1 to reject, start all zeros for callback to work
opt.trl          = trl;
opt.trlop        = 1;
opt.updatethreshold = true;
opt.zmax         = zmax;
opt.zsum         = zsum;

if ~thresholdsum
  opt.zval = zmax;
else
  opt.zval = zsum;
end
opt.zindx = zindx;
if nargin==1
  opt.data = {};
else
  opt.data = data;
end

% if strcmp(cfg.artfctdef.zvalue.interactive, 'yes')
%   set(h, 'visible', 'on');
%   set(h, 'CloseRequestFcn', @cleanup_cb);
%   % give graphical feedback and allow the user to modify the threshold
%   set(h, 'position', [100 200 900 400]);
%   h1 = axes('position', [0.05 0.15 0.4 0.8]);
%   h2 = axes('position', [0.5  0.57  0.45 0.38]);
%   h3 = axes('position', [0.5  0.15  0.45 0.32]);
%   opt.h1           = h1;
%   opt.h2           = h2;
%   opt.h3           = h3;
%   
%   setappdata(h, 'opt', opt);
%   artval_cb(h);
%   redraw_cb(h);
%   
%   % make the user interface elements for the data view
%   uicontrol('tag', 'group1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'stop', 'userdata', 'q')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'downarrow')
%   uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'threshold', 'userdata', 'z')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'uparrow')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'shift+downarrow')
%   uicontrol('tag', 'group1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'artifact','userdata', 'a')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'shift+uparrow')
%   uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'keep trial',   'userdata', 'k')
%   uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'reject part', 'userdata', 'r')
%   uicontrol('tag', 'group3', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'reject full', 'userdata', 'R')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<<', 'userdata', 'shift+leftarrow')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'leftarrow')
%   uicontrol('tag', 'group1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'trial', 'userdata', 't')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'rightarrow')
%   uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>>', 'userdata', 'shift+rightarrow')
%   %uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'userdata', 'ctrl+uparrow')
%   %uicontrol('tag', 'group1', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', 'channel','userdata', 'c')
%   %uicontrol('tag', 'group2', 'parent', h, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'userdata', 'ctrl+downarrow')
%   
%   ft_uilayout(h, 'tag', 'group1', 'width', 0.10, 'height', 0.05);
%   ft_uilayout(h, 'tag', 'group2', 'width', 0.05, 'height', 0.05);
%   ft_uilayout(h, 'tag', 'group3', 'width', 0.12, 'height', 0.05);
%   
%   ft_uilayout(h, 'tag', 'group1', 'style', 'pushbutton', 'callback', @keyboard_cb);
%   ft_uilayout(h, 'tag', 'group2', 'style', 'pushbutton', 'callback', @keyboard_cb);
%   ft_uilayout(h, 'tag', 'group3', 'style', 'pushbutton', 'callback', @keyboard_cb);
%   
%   ft_uilayout(h, 'tag', 'group1', 'retag', 'viewui');
%   ft_uilayout(h, 'tag', 'group2', 'retag', 'viewui');
%   ft_uilayout(h, 'tag', 'group3', 'retag', 'viewui');
%   ft_uilayout(h, 'tag', 'viewui', 'BackgroundColor', [0.8 0.8 0.8], 'hpos', 'auto', 'vpos', 0.005);
%   
%   while opt.quit==0
%     uiwait(h);
%     opt = getappdata(h, 'opt');
%   end
%   
% else

% never interactive
set(h, 'CloseRequestFcn', @cleanup_cb);
% give graphical feedback and allow the user to modify the threshold
set(h, 'position', [100 200 900 400]);
h1 = axes('position', [0.05 0.15 0.4 0.8]);
h2 = axes('position', [0.5  0.57  0.45 0.38]);
h3 = axes('position', [0.5  0.15  0.45 0.32]);
opt.h1           = h1;
opt.h2           = h2;
opt.h3           = h3;

% compute the artifacts given the settings in the cfg
setappdata(h, 'opt', opt);
artval_cb(h);

redraw_cb(h);
[path,filename,~] = fileparts(cfg.saveoutput);
saveas(h,[path filesep filename '.png'], 'png');
  
%end

h   = getparent(h);
opt = getappdata(h, 'opt');

% convert to one long vector
dum = zeros(1,max(opt.trl(:,2)));
for trlop=1:opt.numtrl
  dum(opt.trl(trlop,1):opt.trl(trlop,2)) = opt.artval{trlop};
end
artval = dum;

% find the padded artifacts and put them in a Nx2 trl-like matrix
artbeg = find(diff([0 artval])== 1);
artend = find(diff([artval 0])==-1);
artifact = [artbeg(:) artend(:)];

if strcmp(cfg.artfctdef.zvalue.artfctpeak,'yes')
  cnt=1;
  shift=opt.trl(1,1)-1;
  for tt=1:opt.numtrl
    if tt==1
      tind{tt}=find(artifact(:,2)<opt.trl(tt,2));
    else
      tind{tt}=intersect(find(artifact(:,2)<opt.trl(tt,2)),find(artifact(:,2)>opt.trl(tt-1,2)));
    end
    artbegend=[(artifact(tind{tt},1)-opt.trl(tt,1)+1) (artifact(tind{tt},2)-opt.trl(tt,1)+1)];
    for rr=1:size(artbegend,1)
      [mx,mxnd]=max(opt.zval{tt}(artbegend(rr,1):artbegend(rr,2)));
      peaks(cnt)=artifact(tind{tt}(rr),1)+mxnd-1;
      dssartifact(cnt,1)=max(peaks(cnt)+cfg.artfctdef.zvalue.artfctpeakrange(1)*hdr.Fs,opt.trl(tt,1));
      dssartifact(cnt,2)=min(peaks(cnt)+cfg.artfctdef.zvalue.artfctpeakrange(2)*hdr.Fs,opt.trl(tt,2));
      peaks(cnt)=peaks(cnt)-shift;
      dssartifact(cnt,:)=dssartifact(cnt,:)-shift;
      cnt=cnt+1;
    end
    if tt<opt.numtrl
      shift=shift+opt.trl(tt+1,1)-opt.trl(tt,2)-1;
    end
    clear artbegend
  end
  cfg.artfctdef.zvalue.peaks=peaks';
  cfg.artfctdef.zvalue.dssartifact=dssartifact;
end

% remember the artifacts that were found
cfg.artfctdef.zvalue.artifact = artifact;

% also update the threshold
cfg.artfctdef.zvalue.cutoff   = opt.threshold;

fprintf('detected %d artifacts\n', size(artifact,1));

delete(h);

return
% do the general cleanup and bookkeeping at the end of the function
%ft_postamble provenance
%ft_postamble previous data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function artval_cb(h, eventdata)

opt = getappdata(h, 'opt');

artval = cell(1,opt.numtrl);
for trlop=1:opt.numtrl
  if opt.thresholdsum,
    % threshold the accumulated z-values
    artval{trlop} = opt.zsum{trlop}>opt.threshold;
  else
    % threshold the max z-values
    artval{trlop} = opt.zmax{trlop}>opt.threshold;
  end
  % pad the artifacts
  artbeg = find(diff([0 artval{trlop}])== 1);
  artend = find(diff([artval{trlop} 0])==-1);
  artbeg = artbeg - opt.artpadding;
  artend = artend + opt.artpadding;
  artbeg(artbeg<1) = 1;
  artend(artend>length(artval{trlop})) = length(artval{trlop});
  for artlop=1:length(artbeg)
    artval{trlop}(artbeg(artlop):artend(artlop)) = 1;
  end
  opt.trialok(trlop) = isempty(artbeg);
end

for trlop = find(opt.keep==1 & opt.trialok==0)
  % overrule the objective criterion, i.e. keep the trial when the user
  % wants to keep it
  artval{trlop}(:) = 0;
end

for trlop = find(opt.keep<0 & opt.trialok==1)
  % if the user specifies that the trial is not OK
  % reject the whole trial if there is no extra-threshold data,
  % otherwise use the artifact as found by the thresholding
  if opt.thresholdsum && opt.keep(trlop)==-1,
    % threshold the accumulated z-values
    artval{trlop} = opt.zsum{trlop}>opt.threshold;
  elseif opt.keep(trlop)==-1
    % threshold the max z-values
    artval{trlop} = opt.zmax{trlop}>opt.threshold;
  elseif opt.keep(trlop)==-2
    artval{trlop}(:) = 1;
  end
  % pad the artifacts
  artbeg = find(diff([0 artval{trlop}])== 1);
  artend = find(diff([artval{trlop} 0])==-1);
  artbeg = artbeg - opt.artpadding;
  artend = artend + opt.artpadding;
  artbeg(artbeg<1) = 1;
  artend(artend>length(artval{trlop})) = length(artval{trlop});
  if ~isempty(artbeg)
    for artlop=1:length(artbeg)
      artval{trlop}(artbeg(artlop):artend(artlop)) = 1;
    end
  else
    artval{trlop}(:) = 1;
  end
end

for trlop = find(opt.keep==-2 & opt.trialok==0)
  % if the user specifies the whole trial to be rejected define the whole
  % segment to be bad
  artval{trlop}(:) = 1;
  % pad the artifacts
  artbeg = find(diff([0 artval{trlop}])== 1);
  artend = find(diff([artval{trlop} 0])==-1);
  artbeg = artbeg - opt.artpadding;
  artend = artend + opt.artpadding;
  artbeg(artbeg<1) = 1;
  artend(artend>length(artval{trlop})) = length(artval{trlop});
  if ~isempty(artbeg)
    for artlop=1:length(artbeg)
      artval{trlop}(artbeg(artlop):artend(artlop)) = 1;
    end
  else
    artval{trlop}(:) = 1;
  end
end

opt.artval = artval;
setappdata(h, 'opt', opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyboard_cb(h, eventdata)

if isobject(eventdata)
  % this happens for MATLAB2014b and up, see http://bugzilla.fcdonders.nl/show_bug.cgi?id=2857
  % FIXME the keboard does not seem to work at all at the moment, hence the following work around solves it for now
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
elseif isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(eventdata);
end

% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h = getparent(h);
opt = getappdata(h, 'opt');

switch key
  case 'leftarrow' % change trials
    opt.trlop = max(opt.trlop - 1, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'shift+leftarrow'
    opt.trlop = max(opt.trlop - 10, 1); % should not be smaller than 1
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'rightarrow'
    opt.trlop = min(opt.trlop + 1, opt.numtrl); % should not be larger than the number of trials
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'shift+rightarrow'
    opt.trlop = min(opt.trlop + 10, opt.numtrl); % should not be larger than the number of trials
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'uparrow' % change threshold
    opt.threshold = opt.threshold+0.5;
    opt.updatethreshold = true;
    setappdata(h, 'opt', opt);
    artval_cb(h, eventdata);
    redraw_cb(h, eventdata);
    opt = getappdata(h, 'opt'); % grab the opt-structure from the handle because it has been adjusted in the callbacks
    opt.updatethreshold = false;
    setappdata(h, 'opt', opt);
  case 'downarrow'
    opt.threshold = opt.threshold-0.5;
    opt.updatethreshold = true;
    setappdata(h, 'opt', opt);
    artval_cb(h, eventdata);
    redraw_cb(h, eventdata);
    opt = getappdata(h, 'opt'); % grab the opt-structure from the handle because it has been adjusted in the callbacks
    opt.updatethreshold = false;
    setappdata(h, 'opt', opt);
  case 'shift+uparrow' % change artifact
    artfctindx = find(opt.trialok == 0);
    sel        = find(artfctindx>opt.trlop);
    if ~isempty(sel)
      opt.trlop = artfctindx(sel(1));
    end
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'shift+downarrow'
    artfctindx = find(opt.trialok == 0);
    sel        = find(artfctindx<opt.trlop);
    if ~isempty(sel)
      opt.trlop = artfctindx(sel(end));
    end
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'ctrl+uparrow' % change channel
    if strcmp(opt.channel, 'artifact')
      [dum, indx] = max(opt.zval);
      sgnind      = opt.zindx(indx);
    else
      if ~isempty(opt.data)
        sgnind  = match_str(opt.channel, opt.data.label);
        selchan = match_str(opt.artcfg.channel, opt.channel);
      else
        sgnind  = match_str(opt.channel,   opt.hdr.label);
        selchan = match_str(opt.artcfg.channel, opt.channel);
      end
    end
    numchan = numel(opt.artcfg.channel);
    chansel = min(selchan+1, numchan);
    % convert numeric array into cell-array with channel labels
    opt.channel = tmpchan(chansel);
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'ctrl+downarrow'
    tmpchan = [opt.artcfg.channel;{'artifact'}]; % append the 'artifact' channel
    chansel = match_str(tmpchan, opt.channel);
    chansel = max(chansel-1, 1);
    % convert numeric array into cell-array with channel labels
    opt.channel = tmpchan(chansel);
    setappdata(h, 'opt', opt);
    redraw_cb(h, eventdata);
  case 'a'
    % select the artifact to display
    response = inputdlg(sprintf('artifact trial to display'), 'specify', 1, {num2str(opt.trlop)});
    if ~isempty(response)
      artfctindx = find(opt.trialok == 0);
      sel        = str2double(response);
      sel        = min(numel(artfctindx), sel);
      sel        = max(1,                 sel);
      opt.trlop  = artfctindx(sel);
      setappdata(h, 'opt', opt);
      redraw_cb(h, eventdata);
    end
  case 'c'
    % select channels
    %     select = match_str([opt.artcfg.channel;{'artifact'}], opt.channel);
    %     opt.channel = select_channel_list([opt.artcfg.channel;{'artifact'}], select);
    %     setappdata(h, 'opt', opt);
    %     redraw_cb(h, eventdata);
  case 'q'
    setappdata(h, 'opt', opt);
    cleanup_cb(h);
  case 't'
    % select the trial to display
    response = inputdlg(sprintf('trial to display'), 'specify', 1, {num2str(opt.trlop)});
    if ~isempty(response)
      opt.trlop = str2double(response);
      opt.trlop = min(opt.trlop, opt.numtrl); % should not be larger than the number of trials
      opt.trlop = max(opt.trlop, 1); % should not be smaller than 1
      setappdata(h, 'opt', opt);
      redraw_cb(h, eventdata);
    end
  case 'z'
    % select the threshold
    response = inputdlg('z-threshold', 'specify', 1, {num2str(opt.threshold)});
    if ~isempty(response)
      opt.threshold = str2double(response);
      opt.updatethreshold = true;
      setappdata(h, 'opt', opt);
      artval_cb(h, eventdata);
      redraw_cb(h, eventdata);
      opt = getappdata(h, 'opt'); % grab the opt-structure from the handle because it has been adjusted in the callbacks
      opt.updatethreshold = false;
      setappdata(h, 'opt', opt);
    end
  case 'k'
    opt.keep(opt.trlop) = 1;
    setappdata(h, 'opt', opt);
    artval_cb(h);
    redraw_cb(h);
    opt = getappdata(h, 'opt');
  case 'r'
    % only of the trial contains a partial artifact
    if opt.trialok(opt.trlop) == 0
      opt.keep(opt.trlop) = -1;
    end
    setappdata(h, 'opt', opt);
    artval_cb(h);
    redraw_cb(h);
    opt = getappdata(h, 'opt');
  case 'R'
    opt.keep(opt.trlop) = -2;
    setappdata(h, 'opt', opt);
    artval_cb(h);
    redraw_cb(h);
    opt = getappdata(h, 'opt');
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    setappdata(h, 'opt', opt);
    help_cb(h);
end
uiresume(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function redraw_cb(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

% make a local copy of the relevant variables
trlop     = opt.trlop;
artval    = opt.artval{trlop};
zindx     = opt.zindx{trlop};
zval      = opt.zval{trlop};
cfg       = opt.cfg;
artcfg    = opt.artcfg;
hdr       = opt.hdr;
trl       = opt.trl;
trlpadsmp = round(artcfg.trlpadding*hdr.Fs);
channel   = opt.channel;

% determine the channel with the highest z-value to be displayed
% this is default behaviour but can be overruled in the gui
if strcmp(channel, 'artifact')
  [dum, indx] = max(zval);
  sgnind      = zindx(indx);
else
  if ~isempty(opt.data)
    sgnind = match_str(channel, opt.data.label);
  else
    sgnind = match_str(channel, hdr.label);
  end
end

if ~isempty(opt.data)
  data = ft_fetch_data(opt.data, 'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'));
else
  data = ft_read_data(cfg.datafile,   'header', hdr, 'begsample', trl(trlop,1), 'endsample', trl(trlop,2), 'chanindx', sgnind, 'checkboundary', strcmp(cfg.continuous,'no'));
end
%data = preproc(data, '', hdr.Fs, artcfg, [], artcfg.fltpadding, artcfg.fltpadding);
str  = sprintf('trial %3d, channel %s', opt.trlop, hdr.label{sgnind});
fprintf('showing %s\n', str);

%-----------------------------
% plot summary in left subplot
subplot(opt.h1);hold on;

% plot as a blue line only once
if isempty(get(opt.h1, 'children'))
  for k = 1:opt.numtrl
    xval = opt.trl(k,1):opt.trl(k,2);
    if opt.thresholdsum,
      yval = opt.zsum{k};
    else
      yval = opt.zmax{k};
    end
    plot(opt.h1, xval, yval, 'linestyle', '-', 'color', 'b', 'displayname', 'data');
    xlabel('samples');
    ylabel('z-value');
  end
end
h1children = get(opt.h1, 'children');

% plot trial box
boxhandle = findall(h1children, 'displayname', 'highlight');
if isempty(boxhandle)
  % draw it
  xval = trl(opt.trlop,1):trl(opt.trlop,2);
  if opt.thresholdsum,
    yval = opt.zsum{opt.trlop};
  else
    yval = opt.zmax{opt.trlop};
  end
  plot(opt.h1, xval, yval, 'linestyle', '-', 'color', 'm', 'linewidth', 2, 'displayname', 'highlight');
else
  % update it
  xval = trl(opt.trlop,1):trl(opt.trlop,2);
  if opt.thresholdsum,
    yval = opt.zsum{opt.trlop};
  else
    yval = opt.zmax{opt.trlop};
  end
  set(boxhandle,  'XData', xval);
  set(boxhandle,  'YData', yval);
end

% plot as red lines the suprathreshold data points
thrhandle = findall(h1children, 'displayname', 'reddata');
if isempty(thrhandle)
  % they have to be drawn
  for k = 1:opt.numtrl
    xval = trl(k,1):trl(k,2);
    if opt.thresholdsum,
      yval = opt.zsum{k};
    else
      yval = opt.zmax{k};
    end
    dum = yval<=opt.threshold;
    yval(dum) = nan;
    plot(opt.h1, xval, yval, 'linestyle', '-', 'color', [1 0 0], 'displayname', 'reddata');
  end
  hline(opt.threshold, 'color', 'r', 'linestyle', ':', 'displayname', 'threshline');
elseif ~isempty(thrhandle) && opt.updatethreshold
  % they can be updated
  for k = 1:opt.numtrl
    xval = trl(k,1):trl(k,2);
    if opt.thresholdsum,
      yval = opt.zsum{k};
    else
      yval = opt.zmax{k};
    end
    dum = yval<=opt.threshold;
    yval(dum) = nan;
    set(thrhandle(k), 'XData', xval);
    set(thrhandle(k), 'YData', yval);
  end
  set(findall(h1children, 'displayname', 'threshline'), 'YData', [1 1].*opt.threshold);
end

%--------------------------------------------------
% get trial specific x-axis values and padding info
xval = ((trl(opt.trlop,1):trl(opt.trlop,2))-trl(opt.trlop,1)+trl(opt.trlop,3))./opt.hdr.Fs;
if trlpadsmp>0
  sel    = trlpadsmp:(size(data,2)-trlpadsmp);
  selpad = 1:size(data,2);
else
  sel    = 1:size(data,2);
  selpad = sel;
end

% plot data of most aberrant channel in upper subplot
subplot(opt.h2); hold on
if isempty(get(opt.h2, 'children'))
  % do the plotting
  plot(xval(selpad), data(selpad), 'color', [0.5 0.5 1], 'displayname', 'line1');
  plot(xval(sel),    data(sel),    'color', [0 0 1],     'displayname', 'line2');
  vline(xval(  1)+(trlpadsmp-1/opt.hdr.Fs),     'color', [0 0 0],     'displayname', 'vline1');
  vline(xval(end)-(trlpadsmp/opt.hdr.Fs),       'color', [0 0 0],     'displayname', 'vline2');
  data(~artval) = nan;
  plot(xval, data, 'r-', 'displayname', 'line3');
  xlabel('time(s)');
  ylabel('uV or Tesla');
  xlim([xval(1) xval(end)]);
  title(str);
else
  % update in the existing handles
  h2children = get(opt.h2, 'children');
  set(findall(h2children, 'displayname', 'vline1'), 'visible', 'off');
  set(findall(h2children, 'displayname', 'vline2'), 'visible', 'off');
  set(findall(h2children, 'displayname', 'line1'), 'XData', xval(selpad));
  set(findall(h2children, 'displayname', 'line1'), 'YData', data(selpad));
  set(findall(h2children, 'displayname', 'line2'), 'XData', xval(sel));
  set(findall(h2children, 'displayname', 'line2'), 'YData', data(sel));
  data(~artval) = nan;
  set(findall(h2children, 'displayname', 'line3'),  'XData', xval);
  set(findall(h2children, 'displayname', 'line3'),  'YData', data);
  abc2 = axis(opt.h2);
  set(findall(h2children, 'displayname', 'vline1'), 'XData', [1 1]*xval(  1)+(trlpadsmp-1/opt.hdr.Fs));
  set(findall(h2children, 'displayname', 'vline1'), 'YData', abc2(3:4));
  set(findall(h2children, 'displayname', 'vline2'), 'XData', [1 1]*xval(end)-(trlpadsmp/opt.hdr.Fs));
  set(findall(h2children, 'displayname', 'vline2'), 'YData', abc2(3:4));
  set(findall(h2children, 'displayname', 'vline1'), 'visible', 'on');
  set(findall(h2children, 'displayname', 'vline2'), 'visible', 'on');
  str = sprintf('trial %3d, channel %s', opt.trlop, hdr.label{sgnind});
  title(str);
  xlim([xval(1) xval(end)]);
end

% plot z-values in lower subplot
subplot(opt.h3); hold on;
if isempty(get(opt.h3, 'children'))
  % do the plotting
  plot(xval(selpad), zval(selpad), 'color', [0.5 0.5 1], 'displayname', 'line1b');
  plot(xval(sel),    zval(sel),    'color', [0 0 1],     'displayname', 'line2b');
  hline(opt.threshold, 'color', 'r', 'linestyle', ':', 'displayname', 'threshline');
  vline(xval(  1)+(trlpadsmp-1/opt.hdr.Fs),     'color', [0 0 0],     'displayname', 'vline1b');
  vline(xval(end)-(trlpadsmp/opt.hdr.Fs),       'color', [0 0 0],     'displayname', 'vline2b');
  zval(~artval) = nan;
  plot(xval, zval, 'r-', 'displayname', 'line3b');
  xlabel('time(s)');
  ylabel('z-value');
  xlim([xval(1) xval(end)]);
else
  % update in the existing handles
  h3children = get(opt.h3, 'children');
  set(findall(h3children, 'displayname', 'vline1b'), 'visible', 'off');
  set(findall(h3children, 'displayname', 'vline2b'), 'visible', 'off');
  set(findall(h3children, 'displayname', 'line1b'), 'XData', xval(selpad));
  set(findall(h3children, 'displayname', 'line1b'), 'YData', zval(selpad));
  set(findall(h3children, 'displayname', 'line2b'), 'XData', xval(sel));
  set(findall(h3children, 'displayname', 'line2b'), 'YData', zval(sel));
  zval(~artval) = nan;
  set(findall(h3children, 'displayname', 'line3b'),  'XData', xval);
  set(findall(h3children, 'displayname', 'line3b'),  'YData', zval);
  set(findall(h3children, 'displayname', 'threshline'), 'YData', [1 1].*opt.threshold);
  set(findall(h3children, 'displayname', 'threshline'), 'XData', xval([1 end]));
  abc = axis(opt.h3);
  set(findall(h3children, 'displayname', 'vline1b'), 'XData', [1 1]*xval(  1)+(trlpadsmp-1/opt.hdr.Fs));
  set(findall(h3children, 'displayname', 'vline1b'), 'YData', abc(3:4));
  set(findall(h3children, 'displayname', 'vline2b'), 'XData', [1 1]*xval(end)-(trlpadsmp/opt.hdr.Fs));
  set(findall(h3children, 'displayname', 'vline2b'), 'YData', abc(3:4));
  set(findall(h3children, 'displayname', 'vline1b'), 'visible', 'on');
  set(findall(h3children, 'displayname', 'vline2b'), 'visible', 'on');
  xlim([xval(1) xval(end)]);
end

setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup_cb(h, eventdata)
opt = getappdata(h, 'opt');
opt.quit = true;
setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;

% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
  shiftInd = match_str(eventdata.Modifier, 'shift');
  if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
    % now we now it was a numpad keystroke (numeric character sent AND
    % shift modifier present)
    key = eventdata.Character;
    eventdata.Modifier(shiftInd) = []; % strip the shift modifier
  end
elseif ispc()
  if strfind(eventdata.Key, 'numpad')
    key = eventdata.Character;
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end

function time = offset2time(offset, fsample, nsamples)
% OFFSET2TIME converts the offset of a trial definition into a time-axis
% according to the definition from DEFINETRIAL
%
% Use as
%   [time] = offset2time(offset, fsample, nsamples)
%
offset   = double(offset);
nsamples = double(nsamples);

time = (offset + (0:(nsamples-1)))/fsample;

function [dat, label, time, cfg] = preproc(dat, label, time, cfg, begpadding, endpadding)

% PREPROC applies various preprocessing steps on a piece of EEG/MEG data
% that already has been read from a data file.
%
% compute fsample
fsample = 1./nanmean(diff(time));

if nargin<5 || isempty(begpadding)
  begpadding = 0;
end
if nargin<6 || isempty(endpadding)
  endpadding = 0;
end

if iscell(cfg)
  % recurse over the subsequent preprocessing stages
  if begpadding>0 || endpadding>0
    error('multiple preprocessing stages are not supported in combination with filter padding');
  end
  for i=1:length(cfg)
    tmpcfg = cfg{i};
    if nargout==1
      [dat                     ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==2
      [dat, label              ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==3
      [dat, label, time        ] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
    elseif nargout==4
      [dat, label, time, tmpcfg] = preproc(dat, label, time, tmpcfg, begpadding, endpadding);
      cfg{i} = tmpcfg;
    end
  end
  % ready with recursing over the subsequent preprocessing stages
  return
end

% set the defaults for the rereferencing options
if ~isfield(cfg, 'reref'),        cfg.reref = 'no';             end
if ~isfield(cfg, 'refchannel'),   cfg.refchannel = {};          end
if ~isfield(cfg, 'implicitref'),  cfg.implicitref = [];         end
% set the defaults for the signal processing options
if ~isfield(cfg, 'polyremoval'),  cfg.polyremoval = 'no';       end
if ~isfield(cfg, 'polyorder'),    cfg.polyorder = 2;            end
if ~isfield(cfg, 'detrend'),      cfg.detrend = 'no';           end
if ~isfield(cfg, 'demean'),       cfg.demean  = 'no';           end
if ~isfield(cfg, 'baselinewindow'), cfg.baselinewindow = 'all'; end
if ~isfield(cfg, 'dftfilter'),    cfg.dftfilter = 'no';         end
if ~isfield(cfg, 'lpfilter'),     cfg.lpfilter = 'no';          end
if ~isfield(cfg, 'hpfilter'),     cfg.hpfilter = 'no';          end
if ~isfield(cfg, 'bpfilter'),     cfg.bpfilter = 'no';          end
if ~isfield(cfg, 'bsfilter'),     cfg.bsfilter = 'no';          end
if ~isfield(cfg, 'lpfiltord'),    cfg.lpfiltord = [];           end
if ~isfield(cfg, 'hpfiltord'),    cfg.hpfiltord = [];           end
if ~isfield(cfg, 'bpfiltord'),    cfg.bpfiltord = [];           end
if ~isfield(cfg, 'bsfiltord'),    cfg.bsfiltord = [];           end
if ~isfield(cfg, 'lpfilttype'),   cfg.lpfilttype = 'but';       end
if ~isfield(cfg, 'hpfilttype'),   cfg.hpfilttype = 'but';       end
if ~isfield(cfg, 'bpfilttype'),   cfg.bpfilttype = 'but';       end
if ~isfield(cfg, 'bsfilttype'),   cfg.bsfilttype = 'but';       end
if ~isfield(cfg, 'lpfiltdir'),    if strcmp(cfg.lpfilttype, 'firws'), cfg.lpfiltdir = 'onepass-zerophase'; else cfg.lpfiltdir = 'twopass'; end, end
if ~isfield(cfg, 'hpfiltdir'),    if strcmp(cfg.hpfilttype, 'firws'), cfg.hpfiltdir = 'onepass-zerophase'; else cfg.hpfiltdir = 'twopass'; end, end
if ~isfield(cfg, 'bpfiltdir'),    if strcmp(cfg.bpfilttype, 'firws'), cfg.bpfiltdir = 'onepass-zerophase'; else cfg.bpfiltdir = 'twopass'; end, end
if ~isfield(cfg, 'bsfiltdir'),    if strcmp(cfg.bsfilttype, 'firws'), cfg.bsfiltdir = 'onepass-zerophase'; else cfg.bsfiltdir = 'twopass'; end, end
if ~isfield(cfg, 'lpinstabilityfix'),    cfg.lpinstabilityfix = 'no';    end
if ~isfield(cfg, 'hpinstabilityfix'),    cfg.hpinstabilityfix = 'no';    end
if ~isfield(cfg, 'bpinstabilityfix'),    cfg.bpinstabilityfix = 'no';    end
if ~isfield(cfg, 'bsinstabilityfix'),    cfg.bsinstabilityfix = 'no';    end
if ~isfield(cfg, 'lpfiltdf'),     cfg.lpfiltdf = [];            end
if ~isfield(cfg, 'hpfiltdf'),     cfg.hpfiltdf = [];            end
if ~isfield(cfg, 'bpfiltdf'),     cfg.bpfiltdf = [];            end
if ~isfield(cfg, 'bsfiltdf'),     cfg.bsfiltdf = [];            end
if ~isfield(cfg, 'lpfiltwintype'),cfg.lpfiltwintype = 'hamming';end
if ~isfield(cfg, 'hpfiltwintype'),cfg.hpfiltwintype = 'hamming';end
if ~isfield(cfg, 'bpfiltwintype'),cfg.bpfiltwintype = 'hamming';end
if ~isfield(cfg, 'bsfiltwintype'),cfg.bsfiltwintype = 'hamming';end
if ~isfield(cfg, 'lpfiltdev'),    cfg.lpfiltdev = [];           end
if ~isfield(cfg, 'hpfiltdev'),    cfg.hpfiltdev = [];           end
if ~isfield(cfg, 'bpfiltdev'),    cfg.bpfiltdev = [];           end
if ~isfield(cfg, 'bsfiltdev'),    cfg.bsfiltdev = [];           end
if ~isfield(cfg, 'plotfiltresp'), cfg.plotfiltresp = 'no';      end
if ~isfield(cfg, 'usefftfilt'),   cfg.usefftfilt = 'no';        end
if ~isfield(cfg, 'medianfilter'), cfg.medianfilter  = 'no';     end
if ~isfield(cfg, 'medianfiltord'),cfg.medianfiltord = 9;        end
if ~isfield(cfg, 'dftfreq'),      cfg.dftfreq = [50 100 150];   end
if ~isfield(cfg, 'hilbert'),      cfg.hilbert = 'no';           end
if ~isfield(cfg, 'derivative'),   cfg.derivative = 'no';        end
if ~isfield(cfg, 'rectify'),      cfg.rectify = 'no';           end
if ~isfield(cfg, 'boxcar'),       cfg.boxcar = 'no';            end
if ~isfield(cfg, 'absdiff'),      cfg.absdiff = 'no';           end
if ~isfield(cfg, 'precision'),    cfg.precision = [];           end
if ~isfield(cfg, 'conv'),         cfg.conv = 'no';              end
if ~isfield(cfg, 'montage'),      cfg.montage = 'no';           end
if ~isfield(cfg, 'dftinvert'),    cfg.dftinvert = 'no';         end
if ~isfield(cfg, 'standardize'),  cfg.standardize = 'no';       end
if ~isfield(cfg, 'denoise'),      cfg.denoise = '';             end
if ~isfield(cfg, 'subspace'),     cfg.subspace = [];            end
if ~isfield(cfg, 'custom'),       cfg.custom = '';              end
if ~isfield(cfg, 'resample'),     cfg.resample = '';            end

% test whether the MATLAB signal processing toolbox is available
if strcmp(cfg.medianfilter, 'yes') && ~ft_hastoolbox('signal')
  error('median filtering requires the MATLAB signal processing toolbox');
end

% do a sanity check on the filter configuration
if strcmp(cfg.bpfilter, 'yes') && ...
    (strcmp(cfg.hpfilter, 'yes') || strcmp(cfg.lpfilter,'yes')),
  error('you should not apply both a bandpass AND a lowpass/highpass filter');
end

% do a sanity check on the hilbert transform configuration
if strcmp(cfg.hilbert, 'yes') && ~strcmp(cfg.bpfilter, 'yes')
  error('hilbert transform should be applied in conjunction with bandpass filter')
end

% do a sanity check on hilbert and rectification
if strcmp(cfg.hilbert, 'yes') && strcmp(cfg.rectify, 'yes')
  error('hilbert transform and rectification should not be applied both')
end

% do a sanity check on the rereferencing/montage
if ~strcmp(cfg.reref, 'no') && ~strcmp(cfg.montage, 'no')
  error('cfg.reref and cfg.montage are mutually exclusive')
end

% lnfilter is no longer used
if isfield(cfg, 'lnfilter') && strcmp(cfg.lnfilter, 'yes')
  error('line noise filtering using the option cfg.lnfilter is not supported any more, use cfg.bsfilter instead')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the rereferencing in case of EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.implicitref) && ~any(match_str(cfg.implicitref,label))
  label = {label{:} cfg.implicitref};
  dat(end+1,:) = 0;
end

if strcmp(cfg.reref, 'yes'),
  cfg.refchannel = ft_channelselection(cfg.refchannel, label);
  refindx = match_str(label, cfg.refchannel);
  if isempty(refindx)
    error('reference channel was not found')
  end
  dat = ft_preproc_rereference(dat, refindx);
end

if ~strcmp(cfg.montage, 'no') && ~isempty(cfg.montage)
  % this is an alternative approach for rereferencing, with arbitrary complex linear combinations of channels
  tmp.trial = {dat};
  tmp.label = label;
  tmp = ft_apply_montage(tmp, cfg.montage, 'feedback', 'none');
  dat = tmp.trial{1};
  label = tmp.label;
  clear tmp
end

if any(any(isnan(dat)))
  % filtering is not possible for at least a selection of the data
  warning_once('data contains NaNs, no filtering applied');
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the filtering on the padded data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.denoise),
  hflag    = isfield(cfg.denoise, 'hilbert') && strcmp(cfg.denoise.hilbert, 'yes');
  datlabel = match_str(label, cfg.denoise.channel);
  reflabel = match_str(label, cfg.denoise.refchannel);
  tmpdat   = ft_preproc_denoise(dat(datlabel,:), dat(reflabel,:), hflag);
  dat(datlabel,:) = tmpdat;
end

% The filtering should in principle be done prior to the demeaning to
% ensure that the resulting mean over the baseline window will be
% guaranteed to be zero (even if there are filter artifacts). 
% However, the filtering benefits from the data being pulled towards zero,
% causing less edge artifacts. That is why we start by removing the slow
% drift, then filter, and then repeat the demean/detrend/polyremove.
if strcmp(cfg.polyremoval, 'yes')
  nsamples  = size(dat,2);
  begsample = 1        + begpadding;
  endsample = nsamples - endpadding;
  dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample); % this will also demean and detrend
elseif strcmp(cfg.detrend, 'yes')
  nsamples  = size(dat,2);
  begsample = 1        + begpadding;
  endsample = nsamples - endpadding;
  dat = ft_preproc_polyremoval(dat, 1, begsample, endsample); % this will also demean
elseif strcmp(cfg.demean, 'yes')
  nsamples  = size(dat,2);
  begsample = 1        + begpadding;
  endsample = nsamples - endpadding;
  dat = ft_preproc_polyremoval(dat, 0, begsample, endsample);
end

if strcmp(cfg.medianfilter, 'yes'), dat = ft_preproc_medianfilter(dat, cfg.medianfiltord); end
if strcmp(cfg.lpfilter, 'yes'),     dat = ft_preproc_lowpassfilter(dat, fsample, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir, cfg.lpinstabilityfix, cfg.lpfiltdf, cfg.lpfiltwintype, cfg.lpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
if strcmp(cfg.hpfilter, 'yes'),     dat = ft_preproc_highpassfilter(dat, fsample, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir, cfg.hpinstabilityfix, cfg.hpfiltdf, cfg.hpfiltwintype, cfg.hpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
if strcmp(cfg.bpfilter, 'yes'),     dat = ft_preproc_bandpassfilter(dat, fsample, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir, cfg.bpinstabilityfix, cfg.bpfiltdf, cfg.bpfiltwintype, cfg.bpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
if strcmp(cfg.bsfilter, 'yes')
  for i=1:size(cfg.bsfreq,1)
    % apply a bandstop filter for each of the specified bands, i.e. cfg.bsfreq should be Nx2
    dat = ft_preproc_bandstopfilter(dat, fsample, cfg.bsfreq(i,:), cfg.bsfiltord, cfg.bsfilttype, cfg.bsfiltdir, cfg.bsinstabilityfix, cfg.bsfiltdf, cfg.bsfiltwintype, cfg.bsfiltdev, cfg.plotfiltresp, cfg.usefftfilt);
  end
end
if strcmp(cfg.polyremoval, 'yes')
  % the begin and endsample of the polyremoval period correspond to the complete data minus padding
  nsamples  = size(dat,2);
  begsample = 1        + begpadding;
  endsample = nsamples - endpadding;
  dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample);
end
if strcmp(cfg.detrend, 'yes')
  % the begin and endsample of the detrend period correspond to the complete data minus padding
  nsamples  = size(dat,2);
  begsample = 1        + begpadding;
  endsample = nsamples - endpadding;
  dat = ft_preproc_detrend(dat, begsample, endsample);
end
if strcmp(cfg.demean, 'yes')
  if ischar(cfg.baselinewindow) && strcmp(cfg.baselinewindow, 'all')
    % the begin and endsample of the baseline period correspond to the complete data minus padding
    nsamples  = size(dat,2);
    begsample = 1        + begpadding;
    endsample = nsamples - endpadding;
    dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
  else
    % determine the begin and endsample of the baseline period and baseline correct for it
    begsample = nearest(time, cfg.baselinewindow(1));
    endsample = nearest(time, cfg.baselinewindow(2));
    dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
  end
end
if strcmp(cfg.dftfilter, 'yes')
  datorig = dat;
  for i=1:length(cfg.dftfreq)
    % filter out the 50Hz noise, optionally also the 100 and 150 Hz harmonics
    dat = ft_preproc_dftfilter(dat, fsample, cfg.dftfreq(i));
  end
  if strcmp(cfg.dftinvert, 'yes'),
    dat = datorig - dat;
  end
end
if ~strcmp(cfg.hilbert, 'no')
  dat = ft_preproc_hilbert(dat, cfg.hilbert);
end
if strcmp(cfg.rectify, 'yes'),
  dat = ft_preproc_rectify(dat);
end
if isnumeric(cfg.boxcar)
  numsmp = round(cfg.boxcar*fsample);
  if ~rem(numsmp,2)
    % the kernel should have an odd number of samples
    numsmp = numsmp+1;
  end
  % kernel = ones(1,numsmp) ./ numsmp;
  % dat    = convn(dat, kernel, 'same');
  dat = ft_preproc_smooth(dat, numsmp); % better edge behaviour
end
if isnumeric(cfg.conv)
  kernel = (cfg.conv(:)'./sum(cfg.conv));
  if ~rem(length(kernel),2)
    kernel = [kernel 0];
  end
  dat = convn(dat, kernel, 'same');
end
if strcmp(cfg.derivative, 'yes'),
  dat = ft_preproc_derivative(dat, 1);
end
if strcmp(cfg.absdiff, 'yes'),
  % this implements abs(diff(data), which is required for jump detection
  dat = abs([diff(dat, 1, 2) zeros(size(dat,1),1)]);
end
if strcmp(cfg.standardize, 'yes'),
  dat = ft_preproc_standardize(dat, 1, size(dat,2));
end
if ~isempty(cfg.subspace),
  dat = ft_preproc_subspace(dat, cfg.subspace);
end
if ~isempty(cfg.custom),
  if ~isfield(cfg.custom, 'nargout')
    cfg.custom.nargout = 1;
  end
  if cfg.custom.nargout==1
    dat = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
  elseif cfg.custom.nargout==2
    [dat, time] = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
  end
end
if strcmp(cfg.resample, 'yes')
  if ~isfield(cfg, 'resamplefs')
    cfg.resamplefs = fsample./2;
  end
  if ~isfield(cfg, 'resamplemethod')
    cfg.resamplemethod = 'resample';
  end
  [dat               ] = ft_preproc_resample(dat,  fsample, cfg.resamplefs, cfg.resamplemethod); 
  [time, dum, fsample] = ft_preproc_resample(time, fsample, cfg.resamplefs, cfg.resamplemethod);
end
if ~isempty(cfg.precision)
  % convert the data to another numeric precision, i.e. double, single or int32
  dat = cast(dat, cfg.precision);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the filter padding and do the preprocessing on the remaining trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if begpadding~=0 || endpadding~=0
  dat = ft_preproc_padding(dat, 'remove', begpadding, endpadding);
  if strcmp(cfg.demean, 'yes') || nargout>2
    time = ft_preproc_padding(time, 'remove', begpadding, endpadding);
  end
end

function hline(y, varargin)

% HLINE plot a horizontal line in the current graph

abc = axis;
y = [y y];
x = abc([1 2]);
if length(varargin)==1
  varargin = {'color', varargin{1}};
end
h = line(x, y);
set(h, varargin{:});

function vline(x, varargin)

% VLINE plot a vertical line in the current graph

abc = axis;
x = [x x];
y = abc([3 4]);
if length(varargin)==1
  varargin = {'color', varargin{1}};
end
h = line(x, y);
set(h, varargin{:});



