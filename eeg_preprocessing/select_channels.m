function [selchannelindex, selchannellabels] = select_channels(channellabels,bundlename_or_bundlelabels)
% function [selchannelindex, selchannellabels] = select_channels(channellabels,bundlename_or_bundlelabels)
% select_channels either takes a cell array of electrodes to select, or assumes a standard
% 64-channel BioSemi system for specific bundles. Modify this file to include other systems, or use
% 'ALL_NOSELECTION' to refrain from selecting anything (in this case you have to make a selection
% manually in your EEG files prior to running your analyses if you want to remove EOG etc). This
% function selects channels (e.g. removes HEOG, VEOG or makes some other selection) in a BioSemi and
% other standard 10-20 / 10-05 64-electrode systems. It is very easy to add create other bundles,
% just modify this function to add more bundlenames, e.g. adding a line like:
% bundlenames.YOURBUNDLENAME = {'your' 'electrode' 'selection'};
%
% usage: 
% [channels, channellabels] = select_channels(channellabels,bundlename_or_bundlelabels)
%
% inputs:
%
%   channellabels                   a cell array of the electrodes to select from (i.e. the labels
%                                   of the electrodes that are in the data)
%   bundlename_or_bundlelabels      'ALL_NOSELECTION' (default) determines which channels to select,
%                                   either a bundle name or a cell array of channels. Examples of
%                                   potential values: 'EEG' or 'ALL' selects all EEG channels (NO
%                                   EOG!) 'EOG' selects all EOG channels 'EEG_EOG' selects all EEG
%                                   and all EOG channels 'OCCIP' only occipital, etc: 'PARIET',
%                                   'FRONTAL', 'TEMPORAL', 'OCCIPARIET' or 'CDA' 
%                                   {'Oz', 'Iz', 'O1', 'O2'}; selects OZ, Iz, O1 and O2. The
%                                   resulting bundle name in the results is 'Oz_Iz_O1_O2'.
%
% outputs:
%
%   channels                        indexnumbers of the selected channels channellabels
%   channellabels                   a cell array of the channel names of the selected channels
%
% Internal function of the ADAM toolbox, by J.J.Fahrenfort, 2017, 2018
%
% See also: ADAM_MVPA_FIRSTLEVEL

if nargin<2
    bundlename_or_bundlelabels = 'EEG'; 
end
if strcmpi('ALL',bundlename_or_bundlelabels) 
    bundlename_or_bundlelabels = 'EEG';
end
if iscell(bundlename_or_bundlelabels)
    disp(['Selecting the following channels: ' cell2csv(bundlename_or_bundlelabels) ]);
else
    disp(['Selecting ' bundlename_or_bundlelabels ' channels...']);
end
if strcmpi('EEG',bundlename_or_bundlelabels)
    wraptext('Selecting EEG channels...\nassuming electrodes are defined using labels from the 10-20 system, only electrodes from this system are selected. Edit the select_channels.m function or specify cfg.channelpool = ''ALL_NOSELECTION''; when you do not want to use electrode selection.');
end
% These are the EEG and EOG channels according to the 10-20 / 10-05 definition:
if exist('1005chanlocdata.mat','file')
    load('1005chanlocdata.mat');
else
    chanlocdata = readlocs(trycapfile,'importmode','native'); % from standard 10-20 system
end
all_labels = {chanlocdata(:).labels};
all_types = {chanlocdata(:).type};
bundlenames.EEG =           all_labels(strcmpi(all_types,'EEG'));
bundlenames.EOG =           {all_labels{strcmpi(all_types,'EOG')} 'HEOG_L' 'HEOG_R' 'EOG3'  'EOG4' }; % expand with your own EOG labels if necessary
% These are idiosyncratic bundle definitions, add more bundles or change them to your liking if you want:
bundlenames.OCCIP =         {'PO7'    'PO3'    'O1'     'Iz'    'Oz'    'POz'   'PO8'   'PO4'   'O2'    'PO9'   'PO10'};
bundlenames.PARIET =        {'P1'     'P3'     'P5'     'P7'    'Pz'    'P2'    'P4'    'P6'    'P8' };
bundlenames.FRONTAL =       {'Fp1'    'AF7'    'AF3'    'Fpz'   'Fp2'   'AF8'   'AF4'   'AFz'   'Fz' };
bundlenames.TEMPORAL =      {'FT7'    'C5'     'T7'     'TP7'   'CP5'   'FT8'   'C6'    'T8'    'TP8'   'CP6'};
bundlenames.OCCIPARIET =    {'P1'     'P3'     'P5'     'P7'    'P9'    'Pz'    'P2'    'P4'    'P6'    'P8'    'P10'   'PO7'    'PO3'  'O1'    'Iz'   'Oz'    'POz'   'PO8'   'PO4'   'O2'  'PO9'  'PO10' };
bundlenames.CDA =           {'P5'     'P6'     'P7'     'P8'    'PO7'   'PO8'   'O1'    'O2'    'PO9'    'PO10'};
bundlenames.N2Pc_SPCN =     {'P7'     'P8'     'PO7'    'PO8'   'O1'    'O2'    'PO9'   'PO10'}; % PO3/PO4? scrap O1/O2?
% bundlenames.YOURBUNDLENAME = {'your' 'electrode' 'selection'};

if isempty(channellabels)
    error('No input channel labels were provided, cannot select channels');
end
if iscell(bundlename_or_bundlelabels)
    bundlelabels = bundlename_or_bundlelabels;
elseif strcmpi(bundlename_or_bundlelabels,'ALL_NOSELECTION')
    bundlelabels = channellabels;
elseif strcmpi(bundlename_or_bundlelabels,'EEG_EOG')
    bundlelabels = [ bundlenames.EEG bundlenames.EOG ];
elseif any(strcmpi(bundlename_or_bundlelabels,fieldnames(bundlenames)))
    bundlelabels = bundlenames.(bundlename_or_bundlelabels);
else
    error(['The channel bundle ''' bundlename_or_bundlelabels ''' does not exist, please modify select_channels.m to include it or request one of the existing bundles.']);
end

% finally select labels
selchannelindex = find(ismember(channellabels,bundlelabels));
selchannellabels = channellabels(selchannelindex);
if isempty(selchannellabels)
    error(['could not find any of the channels: ' cell2csv(bundlelabels) ' in the data.']);
end