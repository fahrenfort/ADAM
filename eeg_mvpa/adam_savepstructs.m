function cfg = adam_savepstructs(cfg,varargin)
% ADAM_SAVEPSTRUCTS generates a text file containing all the significant cluster time windows of
% interest. This can be useful to quickly get an overview with p-values.
%
% Use as:
%   adam_savepstructs(cfg, stats1, stats2, ...);
%
% Inputs: 
%   cfg                       = contains input parameters (see below)
%   stats1, stats2, ...       = contains one or more stats variables computed by
%                               adam_compute_group_MVPA, adam_compute_group_ERP, etc.
% 
% The cfg (configuration) input structure can contain the following:
%
% cfg.outputdir               = '' (default) directory path where the text file should be saved
%                               If the cfg struct is empty or does not contain an output directory
%                               specification the function will open a selection window to specify
%                               where the text file should be saved.
% cfg.filename                = 'pstructs.txt' (default) the name of the text file. If this is left
%                               empty it will use the default name 'pstructs.txt'.

if nargin<2
    disp('cannot save significant time windows without stats input, need at least 2 arguments:');
    help adam_savepstructs;
    return
end

% concatenate stats
stats = concat_stats(varargin{:});

% settting some defaults
outputdir = '';
startdir = '';
filename = '';

% unpack config
v2struct(cfg);

% Main routine, is a folder name specified? If not, pop up selection dialog
if isempty(outputdir)
    if ~isfield(cfg,'outputdir')
        cfg.outputdir = '';
    end
    outputdir = uigetdir(startdir,'select directory to save the pstruct results to');
    if ~ischar(outputdir)
        error('no folder was selected');
    end
end
if ~exist(outputdir,'dir')
    error('the specified folder does not exist');
end
cfg.outputdir = outputdir;

% Get filename
if isempty(filename)
    filename = 'pstructs.txt';
end
[~, filename, ~] = fileparts(filename);
filename = [filename '.txt'];
cfg.filename = filename;

% save pstructs
savepstructs(stats,fullfile(outputdir,filename));