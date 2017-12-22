function [selchannelindex, selchannellabels] = select_channels(channellabels,bundlename_or_bundlelabels)
% Either takes a cell array of electrodes to select, or assumes a standard
% 64-channel BioSemi system for specific bundles. Modify this file to
% include other systems, or use 'ALL_NOSELECTION' to refrain from selecting
% anything (in this case you have to make a selection manually in your EEG
% files prior to running your analyses if you want to remove EOG etc). 
% This function selects channels (e.g. removes HEOG, VEOG or makes some
% other selection) in a BioSemi and other standard 10-20 64-electrode
% systems. It is very easy to add create other bundles, just modify this
% function to add more bundlenames, e.g. adding a line like:
% bundlenames.YOURBUNDLENAME = {'your' 'electrode' 'selection'};
%
% usage: 
% [channels, channellabels] = select_channels(channellabels,bundlename_or_bundlelabels)
%
% inputs:
% channellabels: a cell array of the electrodes to select from (i.e. the labels of the electrodes that are in the data)
% bundlename_or_bundlelabels determines which channels to select, either a bundle name:
% bundlename_or_bundlelabels = 'EEG' or 'ALL' takes all EEG channels (NO EOG!)
% bundlename_or_bundlelabels = 'EOG' all EOG channels
% bundlename_or_bundlelabels = 'EEG_EOG' all EEG and all EOG channels
% bundlename_or_bundlelabels = 'OCCIP' only occipital, 'PARIET' only parietal, etc: 'FRONTAL' 'TEMPORAL' 'OCCIPARIET' 'CDA' 
% or a cell array of channels:
% bundlename_or_bundlelabels = {'Oz', 'Iz', 'O1', 'O2'};
%
% outputs:
% channels are the indexnumbers of the selected channels channellabels
% channellabels is a cell array of the channel names of the selected channels
%
% Part of the ADAM toolbox, J.J.Fahrenfort, VU 2016, 2017
if nargin<2
    bundlename_or_bundlelabels = 'EEG'; 
end
if strcmpi('ALL',bundlename_or_bundlelabels) 
    bundlename_or_bundlelabels = 'EEG';
end
if iscell(bundlename_or_bundlelabels)
    disp(['selecting the following channels: ' cellarray2csvstring(bundlename_or_bundlelabels) ]);
else
    disp(['selecting ' bundlename_or_bundlelabels ' channels...']);
end
if strcmpi('EEG',bundlename_or_bundlelabels)
    wraptext('Removing EOG channels...\n\nassuming this is a 64 channel BioSemi 10-20 system, electrodes outside the 64 channel range are always removed. Edit the select_channels.m function if this is undesirable behavior.');
end
% This is the bundle definition, add more bundles to your liking:
bundlenames.EEG =           {'Fp1'    'AF7'    'AF3'    'F1'    'F3'    'F5'    'F7'    'FT7'   'FC5'   'FC3'   'FC1'   'C1'    'C3'    'C5'    'T7'    'TP7'    'CP5'    'CP3' 'CP1'    'P1'    'P3'    'P5'    'P7'    'P9'    'PO7'    'PO3'    'O1'    'Iz'  'Oz'    'POz'    'Pz'    'CPz'    'Fpz'    'Fp2'    'AF8'    'AF4'    'AFz' 'Fz'    'F2'    'F4'    'F6'    'F8'    'FT8'    'FC6'    'FC4'    'FC2'  'FCz'    'Cz'    'C2'    'C4'    'C6'    'T8'    'TP8'    'CP6'    'CP4'  'CP2'    'P2'    'P4'    'P6'    'P8'    'P10'    'PO8'    'PO4'    'O2'  'PO9'    'PO10'};
bundlenames.EOG =           {'VEOG'   'HEOG'   'HEOG_L' 'HEOG_R' 'EOG1' 'EOG2'  'EOG3'  'EOG4'  'LO1'   'LO2'   'IO1'   'IO2'   'SO1'   'SO2' };
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
    error(['could not find any of the channels: ' cellarray2csvstring(bundlelabels) ' in the data.']);
end