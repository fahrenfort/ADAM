function adam_detrend_and_epoch(cfg)
% function adam_detrend_and_epoch(cfg)
% Reads in CONTINUOUS (!) EEGLAB data, removes faulty electrodes, detrends using polynomials, epochs
% and writes out as EPOCHED (!) EEGLAB file. Note that the data are high-pass filtered early on to
% identify faulty electrodes, but that the highpass filtered data is discarded after identifying
% these electrodes. Polynomial detrending is applied to remove drifts in raw data. Faulty
% electrodes are interpolated back in after detrending. filenames: names of files either in a cell
% array or as comma separated string. Wildcards * and ? can be used, e.g. filenames = '*' will take
% all the .set files in the input filepath as sources. The output of this function can serve as
% input for ADAM_MVPA_FIRSTLEVEL.
%
%       cfg.datadir           = string specifiying the directory where the input files are
%                               located;
%       cfg.outputdir         = string specifying the directory where the results should be
%                               saved. Choose an informative name for the analysis.
%       cfg.filenames         = N by 1 cell array containing N strings, in which each string
%                               contains the base of a filename (e.g. cfg.filenames =
%                               {'subject1', 'subject2'};). Do not add a file extension to the
%                               file name, ADAM will automatically look for EEGLAB .set files.
%       cfg.start_epoch       = in seconds (relative to target events)
%       cfg.end_epoch         = in seconds (relative to target events)
%       cfg.polynomial_order  = order of the polynomial that is used for detrending (Default: 30,
%                               but note that higher orders can improve the data even more).
%       cfg.pad_length        = wide epoch padding window used to fit the polynomial on, e.g.
%                               pad_length of 100 creates 50 second pads around both sides of a
%                               trial (Default: 50). This padding window is only used during
%                               detrending and is discarded after the trial itsel is epoched out
%                               (when doing the final epoching).
%       cfg.start_mask        = in seconds, specifies when your 'cognitive' or other events of
%                               interest start
%       cfg.end_mask          = in seconds, specifies when your 'cognitive' or other events of
%                               interest end
%       cfg.mask_only_current = 'yes' or 'no' (default: 'yes'). Specifies whether the algorithm for
%                               any given trial masks out all other trials from the wide window pad,
%                               or only the currently relevant trial. If you have sufficiently long
%                               intertrial intervals you may also set this to false.
%       cfg.event_codes       = array containing the relevant event values for the experiment.
%       cfg.remove_bad_chans  = identifies and interpolates bad electrodes.
%
% The function also saves a graph for each subject in the output folder, containg (from top to
% bottom):
%       1.      A graph of the middle trial in the experiment, plotting the five electrodes with the
%               largest standard deviation across the wide epoch, their 1st order polynomial fits as
%               dotted lines. The black line at the top shows the masks for these five electrodes
%               (in which the gap in the line is the temporal window corresponding to cognitive
%               events that are masked out)
%       2.      The same five electrodes, but now with their 1st order polynomial fits removed,
%               their nth order polynomials (n is specified by polynomial_order) as dotted lines,
%               plus at the top the masks for these five electrodes (same as above).
%       3.      The same five electrodes with their nth order polynomials removed.
%       4.      The butterfly ERP before polynomial removal (baselined at <0 for illustration
%               purposes).
%       5.      The butterfly ERP after polynomial removal (baselined at <0 for illustration
%               purposes, the output data of this function is not baselined).
%
% Usage example:
%
% cfg = [];
% cfg.datadir           = 'path/to/continous_data;
% cfg.outputdir         = 'path/to/detrended_epoched_data';
% cfg.filenames         = {'subject01', 'subject02', 'subject03', 'subject04', 'subject05'};
% cfg.start_epoch       = -.5;
% cfg.end_epoch         = 6;
% cfg.start_mask        = 0;
% cfg.end_mask          = 5.5;
% cfg.polynomial_order  = 30;
% cfg.pad_length        = 50;
% cfg.mask_only_current = 'no';
% cfg.event_codes       = [1001:1024 1101:1124];
% cfg.remove_bad_chans  = true;
%
% adam_detrend_eeg_and_epoch(cfg);
%
% part of the ADAM toolbox, by J.J.Fahrenfort, VU, 2018, 2019
%
% See also ADAM_MVPA_FIRSTLEVEL, ADAM_COMPUTE_GROUP_MVPA, ADAM_PLOT_MVPA

% default values
polynomial_order = 30;
pad_length = 50;
mask_only_current = true;
remove_bad_chans = true;

% unpack cfg
v2struct(cfg);

% do some checking
if ~exist('datadir','var')
    error('You need to specify in cfg.datadir where the data are located');
end
if ~exist('outputdir','var')
    error('You need to specify in cfg.outputdir where the results should be stored');
end
if ~exist('filenames','var')
    error('You need to specify a cell array in cfg.filenames containing the filenames containing the data of each subject');
end
if ~exist('event_codes','var')
    error('You need to specify in cfg.event_codes which event values you want to epoch on');
end
if ~exist('start_epoch','var')
    error('You need to specify in cfg.start_epoch defining the onset time of your epoch (in seconds)');
end
if ~exist('end_epoch','var')
    error('You need to specify in cfg.end_epoch defining the offset time of your epoch (in seconds)');
end
if ~exist('start_mask','var')
    error('You need to specify in cfg.start_mask defining the onset time of your mask (in seconds)');
end
if ~exist('end_mask','var')
    error('You need to specify in cfg.end_mask defining the offset time of your mask (in seconds)');
end
if ~ischar(mask_only_current)
    if mask_only_current
        mask_only_current = 'yes';
    else
        mask_only_current = 'no';
    end
end
if ~ischar(remove_bad_chans)
    if remove_bad_chans
        remove_bad_chans = 'yes';
    else
        remove_bad_chans = 'no';
    end
end
if ~ischar(event_codes)
    event_codes = cond_string(event_codes);
end

% run analysis
if ~exist('qsub','var') || isempty(qsub) % run local
    for cSubj = 1:numel(filenames)
        detrend_and_epoch(datadir,filenames{cSubj},outputdir, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, remove_bad_chans, event_codes);
    end
else % or create qsub files
    create_slurm_files(qsub.functionpath,'detrend_and_epoch',qsub, datadir, filenames, outputdir, start_epoch, end_epoch, polynomial_order, pad_length, start_mask, end_mask, mask_only_current, remove_bad_chans, event_codes);
end

