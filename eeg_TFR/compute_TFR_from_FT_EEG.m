function [TFR, groupTFR] = compute_TFR_from_FT_EEG(FT_EEG,condSet,resample_to,method,tf_baseline,erp_baseline,frequencies)
% compute_TFR_from_FT_EEG computes time frequency representation (TFR) on FT_EEG and returns result.
% resample_to specifies the sampling rate to which to downsample. If left empty or specified as 0,
% no resampling is done. tf_and_erp_baseline is the period that is used as baseline normalization,
% in seconds e.g.: [-.4 -.15]. If this field is 0 or empty, no baseline correction is applied. If
% separated by a semicolon, the first part is for the tf_baseline and the second part for the
% erp_baseline e.g. tf_and_erp_baseline = '-.35,-.1;-.25,0' uses a baseline of -350 to -100 ms for
% the tf baseline and a baseline of -250 to 0 ms for the erp baseline. To apply only one type of
% baseline, set the other type to 0,0, like this: tf_and_erp_baseline = '-.35,-.1;0,0'; applies only
% a TF baseline. You can specify the kind of TF baseline that is used by specifying method: e.g.
% method = 'absolute' (default). Other options: method = 'relative', 'relchange', 'normchange' or
% 'db' method is also used to specify the method(s) used for computing the TFR: method = 'total'
% (default), 'evoked' or 'induced' total and induced can either be used  with 'bin' or 'single',
% e.g.: method = 'relchange,induced,bin', meaning that trials will either be binned per condition
% (see below) or they will be kept as individual trials. Evoked is always the result of binning,
% otherwise there would be no evoked signal to compute the TFR on. Finally, induced can be computed
% using two methods: either by subtracting the average ERP from the bin class condition from each
% trial in that class (method = 'subtr_bin') or by subtracting the average ERP from the original
% condition label from each trial that has that label (method = 'subtr_indiv', default). TFR
% baseline correction will also be applied according the subtr_indiv method, either correcting the
% TFR baseline per newly defined stimulus class, or per original condition label.
% The following names will be added to the TFR in the field fname, depending on the
% settings used to compute it:
% TOT_SINGLE:               Computed individual total power values (computing power for each trial,
%                           no binning, but recoding condition values to new condition numbers in
%                           case bins were specified in the condition set)
% TOT_BIN:                  Computed binned power values (averaging bins after power power
%                           computation)
% EVK_BIN:                  Computed evoked power values (averaging bins
%                           before power computation)
% INDC_SINGLE_SUBTR_BIN:    Computed induced power per trial (each stimulus class' bin average erp
%                           subtracted from each respective trial before power computation)
% INDC_SINGLE_SUBTR_INDIV:  Computed induced power values (each original condition label average erp
%                           subtracted from each respective trial before power computation)
% INDC_BIN_SUBTR_BIN:       Computed induced power values (each stimulus class' bin average erp
%                           subtracted from each respective trial before power computation, binning
%                           after power computation)
% INDC_BIN_SUBTR_INDIV:     Computed induced power values (each original condition label average erp
%                           subtracted from each respective trial before power computation, binning
%                           after power computation)
% frequencies indicates which frequencies should be computed, can be
% indicated as '1:20' or as '1:2:20' or even as '1,3,5,7'. If left empty,
% all frequencies from 2 to 100 are computed in steps of 2 (so the
% equivalent of frequencies = '2:2:100'
% condSet (varargin) can be set up as follows (example):
% condSet{1} = [ 1, 2, 3];
% condSet{2} = [ 4, 5, 6];
% or as a list of comma separated values:
% cond1 = '1,2,3';
% cond2 = '4,5,6';
% Averages single instances from condition 1,2 and 3 into a 'new' instance of condition 1, and
% averages single instances from condition 4, 5 and 6 into a single new instance of condition 2.
% 'Leftover' trials are discarded. See compute_bins_on_FT_EEG for more details.
%
% Internal function of the ADAM toolbox by J.J.Fahrenfort, VU, 2015, 2016, 2018
%
% See also: ADAM_MVPA_FIRSTLEVEL, CLASSIFY_TFR_FROM_EEGLAB_DATA, COMPUTE_TFR, COMPUTE_BINS_ON_FT_EEG

% bunch of sanity checks
if nargin < 7
    disp('you did not specify all the required arguments. See help for more information.');
    return;
end
groupTFR = []; % intialize
if ischar(resample_to)
    resample_to = string2double(resample_to);
end
if isempty(resample_to) 
    resample_to = 0;
end
% check whether random permutation and/or binning should be performed
methods = regexp(method,',','split');
method = 'total';
binning = false;
subtr_method = 'subtr_indiv';
bl_method = 'relchange';
only_group_data = false;
use_splines = false;
for c=1:numel(methods)
    if sum(strcmpi(methods{c},{'total', 'evoked', 'induced'})) == 1
        method = methods{c};
    end
    if strcmpi(methods{c},'single')
        binning = false;
    end
    if strcmpi(methods{c},'bin')
        binning = true;
    end   
    if sum(strcmpi(methods{c},{'absolute', 'relative', 'relchange', 'vssum', 'db'})) == 1
        bl_method = methods{c};
    end
    if sum(strcmpi(methods{c},{'subtr_bin', 'subtr_indiv'})) == 1
        subtr_method = methods{c};
    end
    if sum(strcmpi(methods{c},{'no_individual_trials', 'only_group'})) == 1
        only_group_data = true;
    end
    if sum(strcmpi(methods{c},{'splines','spline'}))
        use_splines = true;
    end
end
if ischar(frequencies)
    frequencies = str2num(frequencies);
end
if isempty(condSet)
    error('cannot find usable trigger specification');
end

% erp_baseline
cfg = [];
cfg.baseline = erp_baseline;
cfg.channel = 'all';
cfg.parameter = 'trial';
FT_EEG = ft_timelockbaseline(cfg,FT_EEG);
if strcmpi(erp_baseline,'no') && binning == true && (strcmp(method,'induced') || strcmp(method,'evoked'))
    disp('WARNING: you should probably ERP baseline your data before you compute evoked or induced bins on the TF data! re-run and specify an ERP baseline.');
end

% if method is evoked, compute evoked data
if strcmp(method,'evoked')
    % compute 'evoked' part of signal
    FT_EVOKED = compute_bins_on_FT_EEG(FT_EEG,condSet,'trial','original');
    % turn into splines
    if use_splines
        FT_EVOKED = compute_spline_on_FT_EEG(FT_EVOKED);
    end
    [TFR, groupTFR] = compute_TFR(FT_EVOKED,resample_to,frequencies);
    clear FT_EVOKED;
    % TFR baseline
    if strcmpi(tf_baseline,'no')
        disp('no TF baseline applied');
    else
        TFR = baseline_TFR(TFR, condSet, tf_baseline, bl_method, subtr_method);
    end
    fname = 'EVK_BIN';
end

% if method is induced
if strcmp(method,'induced')
    % compute 'evoked' part of signal
    if strcmp(subtr_method,'subtr_bin')
        % FT_EVOKED = compute_bins_on_FT_EEG(FT_EEG,condSet);
        FT_EVOKED = compute_erp_on_FT_EEG(FT_EEG,condSet,'trial','bin');
    elseif strcmp(subtr_method,'subtr_indiv')
         % subtracts the erp from each condition in a condSet
        FT_EVOKED = compute_erp_on_FT_EEG(FT_EEG,condSet,'trial','indiv');
    end
    % turn into splines
    if use_splines
        FT_EVOKED = compute_spline_on_FT_EEG(FT_EVOKED);
    end
    % subtract out the evoked responses or the splines
    FT_INDUCED = subtract_evoked_from_FT_EEG(FT_EEG,FT_EVOKED);
    clear FT_EVOKED;
    % compute the TFR
    [TFR, groupTFR] = compute_TFR(FT_INDUCED,resample_to,frequencies);
    clear FT_INDUCED;
    % TFR baseline
    if strcmpi(tf_baseline,'no')
        disp('no TF baseline applied');
        if binning
            disp('WARNING: you might want to apply a TF baseline before binning');
        end
    else
        TFR = baseline_TFR(TFR, condSet, tf_baseline, bl_method, subtr_method);
    end
    if binning
        % compute bins on induced (after power computation) and save
        TFR = compute_bins_on_FT_EEG(TFR,condSet,'powspctrm','original');
        fname = ['INDC_BIN_' upper(subtr_method)];
    else
        fname = ['INDC_SINGLE_' upper(subtr_method)];
    end
end

% if method is total
if strcmp(method,'total')
    % compute the TFR
    [TFR, groupTFR] = compute_TFR(FT_EEG,resample_to,frequencies);
    % TFR baseline
    if strcmpi(tf_baseline,'no')
        disp('no TF baseline applied');
        if binning
            disp('WARNING: you might want to apply a TF baseline before binning');
        end
    else
        TFR = baseline_TFR(TFR, condSet, tf_baseline, bl_method, subtr_method);
    end
    if binning
        % compute bins on total and save
        TFR = compute_bins_on_FT_EEG(TFR,condSet,'powspctrm','original');
        fname = 'TOT_BIN';
    else
        fname = 'TOT_SINGLE';
    end
end

% return the single subject condition averages for group analysis visualizations
TFR.fname = fname;