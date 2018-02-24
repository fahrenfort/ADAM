function FT_EEG = eeglab2ft(EEG,filepath)
% function FT_EEG = eeglab2ft(EEG,filepath)
% Wrapper function to import EEG lab data to a fieldtrip, including
% conditions
% EEG can either be an eeglab struct, or be the filename which contains the
% eeglab data. In the latter case, filepath must also be specified
%
% J.J.Fahrenfort, VU 2014, 2016

if ~isstruct(EEG)
    if isempty(filepath)
        error('if you want to read in the data, you need to specify the path');
    end
    % remove .set if present
    EEG(strfind(EEG,'.set'):end) = [];
    % read data from EEGlab
    EEG = pop_loadset('filename',[EEG '.set'],'filepath',filepath);
end
FT_EEG = eeglab2fieldtrip_local(EEG, 'preprocessing');
% also get events
for cEvents = 1:numel(EEG.epoch)
    if iscell(EEG.epoch(cEvents).eventlatency)
        index = cell2mat(EEG.epoch(cEvents).eventlatency) == 0;
    else
        index = EEG.epoch(cEvents).eventlatency == 0;
    end
    event = EEG.epoch(cEvents).eventtype(index);
    if iscell(event)
        event = event{1};
    end
    if ischar(event)
        event = string2double(event);
        if numel(event) > 1
            warning('cannot convert event to a single numeric value, taking the first value');
            event = event(1);
        end
        if isempty(event) || isnan(event)
            event = NaN;
            warning('event value is NaN when converted to a numeral');
        end
    end
    FT_EEG.trialinfo(cEvents,1) = event;
end
clear EEG;

function data = eeglab2fieldtrip_local(EEG, fieldbox, transform)

if nargin < 2
    error('missing 2nd argument');
    return;
end;

% start with an empty data object 
data = [];

% add the objects that are common to all fieldboxes
tmpchanlocs  = EEG.chanlocs;
data.label   = { tmpchanlocs(EEG.icachansind).labels };
data.fsample = EEG.srate;

% get the electrode positions from the EEG structure: in principle, the number of 
% channels can be more or less than the number of channel locations, i.e. not 
% every channel has a position, or the potential was not measured on every
% position. This is not supported by EEGLAB, but it is supported by FIELDTRIP.

if strcmpi(fieldbox, 'chanloc_withfid')
    % insert "no data channels" in channel structure
    % ----------------------------------------------
    if isfield(EEG.chaninfo, 'nodatchans') && ~isempty( EEG.chaninfo.nodatchans )
        chanlen = length(EEG.chanlocs);
        fields = fieldnames( EEG.chaninfo.nodatchans );
        for index = 1:length(EEG.chaninfo.nodatchans)
            ind = chanlen+index;
            for f = 1:length( fields )
                EEG.chanlocs = setfield(EEG.chanlocs, { ind }, fields{f}, ...
                                        getfield( EEG.chaninfo.nodatchans, { index },  fields{f}));
            end;
        end;
    end;
end;

data.elec.pnt   = zeros(length( EEG.chanlocs ), 3);
for ind = 1:length( EEG.chanlocs )
    data.elec.label{ind} = EEG.chanlocs(ind).labels;
    if ~isempty(EEG.chanlocs(ind).X)
        data.elec.pnt(ind,1) = EEG.chanlocs(ind).X;
        data.elec.pnt(ind,2) = EEG.chanlocs(ind).Y;
        data.elec.pnt(ind,3) = EEG.chanlocs(ind).Z;
    else
        data.elec.pnt(ind,:) = [0 0 0];
    end;
end;

if nargin > 2
    if strcmpi(transform, 'dipfit') 
        if ~isempty(EEG.dipfit.coord_transform)
            disp('Transforming electrode coordinates to match head model');
            transfmat = traditionaldipfit(EEG.dipfit.coord_transform);
            data.elec.pnt = transfmat * [ data.elec.pnt ones(size(data.elec.pnt,1),1) ]';
            data.elec.pnt = data.elec.pnt(1:3,:)';
        else
            disp('Warning: no transformation of electrode coordinates to match head model');
        end;
    end;
end;
        
switch fieldbox
  case 'preprocessing'
    for index = 1:EEG.trials
      data.trial{index}  = EEG.data(:,:,index);
      data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
    end;
    data.label   = { tmpchanlocs(1:EEG.nbchan).labels };

    
  case 'timelockanalysis'
    data.avg  = mean(EEG.data, 3);   
    data.var  = std(EEG.data, [], 3).^2;   
    data.time = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
    data.label   = { tmpchanlocs(1:EEG.nbchan).labels };
    
  case 'componentanalysis'
      if isempty(EEG.icaact)
          icaacttmp = eeg_getica(EEG);
      end
    for index = 1:EEG.trials
      % the trials correspond to the raw data trials, except that they
      % contain the component activations
      try
          if isempty(EEG.icaact)
              data.trial{index} = icaacttmp(:,:,index); % Using icaacttmp to not change EEG structure
          else
              data.trial{index}  = EEG.icaact(:,:,index);
          end
      catch
          
      end;
      data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
    end;
    data.label = [];
    for comp = 1:size(EEG.icawinv,2)
      % the labels correspond to the component activations that are stored in data.trial
      data.label{comp} = sprintf('ica_%03d', comp);
    end
    % get the spatial distribution and electrode positions
    tmpchanlocs    = EEG.chanlocs;
    data.topolabel = { tmpchanlocs(EEG.icachansind).labels };
    data.topo      = EEG.icawinv;
    
  case { 'chanloc' 'chanloc_withfid' }
 
  case 'freqanalysis'
    error('freqanalysis fieldbox not implemented yet')
    
  otherwise
    error('unsupported fieldbox') 
end

try
  % get the full name of the function
  data.cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  data.cfg.version.name = st(i);
end

% add the version details of this function call to the configuration
data.cfg.version.id   = '$Id: eeglab2fieldtrip_local.m,v 1.6 2009-07-02 23:39:29 arno Exp $';

