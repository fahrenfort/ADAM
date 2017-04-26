function [channels, channelnames] = select_channels(channelnames,selectbundle)
% function [channels, channelnames] = select_channels(channelnames,selectbundle)
% selects channels (e.g. removes HEOG, VEOG or makes some other selection)
% selectbundle determines which channels to select:
% selectbundle = 'EEG' or 'ALL' takes all EEG channels (NO EOG!)
% selectbundle = 'EOG' all EOG channels
% selectbundle = 'EEG_EOG' all EEG and all EOG channels
% selectbundle = 'OCCIP', only occipital channels etc 
% outputs:
% channels are the indexnumbers of the selected channels
% channelnames is a cell array of the channel names of the selected channels
% J.J.Fahrenfort, VU 2016
if nargin<2
    selectbundle = 'EEG'; 
end
if strcmpi('ALL',selectbundle) 
    selectbundle = 'EEG';
end
disp(['selecting ' selectbundle ' channels...']);
if strcmpi('EEG',selectbundle)
    wraptext('Removing EOG channels...\n\nassuming this is a 64 channel BioSemi 10-20 system, electrodes outside the 64 channel range are always removed. Edit the select_channels.m function if this is undesirable behavior.');
end
elecSetNames = {'EEG' 'EOG' 'OCCIP' 'PARIET' 'FRONTAL' 'TEMPORAL' 'OCCIPARIET' 'CDA' 'N2Pc_SPCN'};
elecSets = {
    {'Fp1'    'AF7'    'AF3'    'F1'    'F3'    'F5'    'F7'    'FT7'   'FC5'   'FC3'   'FC1'   'C1'    'C3'    'C5'    'T7'    'TP7'    'CP5'    'CP3' 'CP1'    'P1'    'P3'    'P5'    'P7'    'P9'    'PO7'    'PO3'    'O1'    'Iz'  'Oz'    'POz'    'Pz'    'CPz'    'Fpz'    'Fp2'    'AF8'    'AF4'    'AFz' 'Fz'    'F2'    'F4'    'F6'    'F8'    'FT8'    'FC6'    'FC4'    'FC2'  'FCz'    'Cz'    'C2'    'C4'    'C6'    'T8'    'TP8'    'CP6'    'CP4'  'CP2'    'P2'    'P4'    'P6'    'P8'    'P10'    'PO8'    'PO4'    'O2'  'PO9'    'PO10'}, ...
    {'VEOG'   'HEOG'   'HEOG_L' 'HEOG_R' 'EOG1' 'EOG2'  'EOG3'  'EOG4'  'LO1'   'LO2'   'IO1'   'IO2'   'SO1'   'SO2' }, ...
    {'PO7'    'PO3'    'O1'     'Iz'    'Oz'    'POz'   'PO8'   'PO4'   'O2'    'PO9'   'PO10'}, ...
    {'P1'     'P3'     'P5'     'P7'    'Pz'    'P2'    'P4'    'P6'    'P8' }, ...
    {'Fp1'    'AF7'    'AF3'    'Fpz'   'Fp2'   'AF8'   'AF4'   'AFz'   'Fz' }, ...
    {'FT7'    'C5'     'T7'     'TP7'   'CP5'   'FT8'   'C6'    'T8'    'TP8'   'CP6'}, ...
    {'P1'     'P3'     'P5'     'P7'    'P9'    'Pz'    'P2'    'P4'    'P6'    'P8'    'P10'   'PO7'    'PO3'  'O1'    'Iz'   'Oz'    'POz'   'PO8'   'PO4'   'O2'  'PO9'  'PO10' }, ...
    {'P5'     'P6'     'P7'     'P8'    'PO7'   'PO8'   'O1'    'O2'    'PO9'    'PO10'}, ...
    {'P7'     'P8'     'PO7'    'PO8'   'O1'    'O2'     'PO9'    'PO10'}, ... % PO3/PO4? scrap O1/O2?
};
if isempty(channelnames)
    channelnames = [elecSets{strcmpi('EEG',elecSetNames)}  elecSets{strcmpi('EOG',elecSetNames)}];
end
if strcmpi(selectbundle,'EEG_EOG')
    channel_bundle = [elecSets{strcmpi('EEG',elecSetNames)}  elecSets{strcmpi('EOG',elecSetNames)}];
else
    if ~any(strcmpi(selectbundle,elecSetNames))
        error(['The electrode/channel bundle ''' selectbundle ''' does not exist, please modify select_channels.m to include it or request one of the existing bundles.']);
    else
        channel_bundle = elecSets{strcmpi(selectbundle,elecSetNames)};
    end
end
channels = [];
for cChan = 1:numel(channelnames)
    if any(strcmpi(channelnames{cChan},channel_bundle))
        channels = [channels cChan];
    end
end
channelnames = channelnames(channels);