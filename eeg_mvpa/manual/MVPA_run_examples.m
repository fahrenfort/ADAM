% settings
filepath = '/Users/VU-MBP/data';
filename = 'subj01.set';
outpath = '/Users/VU-MBP/MVPA_RESULTS_RAW/';
nFolds = 10;
cChannelSet = 1;
method = 'linear,hr-far';
crossclass = '1';
erp_baseline = '-.25,0';
cond1 = '1,2,3'; cond2 = '4,5,6'; cond3 = '7,8,9'; 
% go
classify_RAW_eeglab_data(filepath,filename,outpath,nFolds,cChannelSet,...
    method,crossclass,erp_baseline,cond1,cond2,cond3);





attended = 1:10;
non_attended = 11:20;
masked = [1:5 11:15];
unmasked = [6:10 16:20];
cond1 = cond_string(attended, masked);
cond2 = cond_string(attended, unmasked);

% settings
filepath = '/Users/VU-MBP/data';
filename = 'subj01.set';
outpath = '/Users/VU-MBP/MVPA_RESULTS_TFR/';
% TFR settings
resample_to = 64;
TFRmethod = 'total,bin';
tf_and_erp_baseline = '0,0;-.25,0';
frequencies = '2:2:30';
% MVPA settings
nFolds = 10;
channelset = 1;
MVPAmethod = 'linear';
crossclass = '0';
cond1 = '1,2,3'; cond2 = '4,5,6'; cond3 = '7,8,9';
% go
classify_TFR_data_from_eeglab(filepath,filename,outpath,resample_to, ...
    TFRmethod,tf_and_erp_baseline,frequencies,nFolds,channelset, ...
    MVPAmethod,crossclass,cond1,cond2, cond3);


% settings
filepath = '/Users/VU-MBP/data';
filename = 'subj01train,subj01test';
outpath = '/Users/VU-MBP/MVPA_RESULTS_RAW/';
nFolds = 1;
cChannelSet = 1;
method = 'linear';
crossclass = '1';
erp_baseline = '-.25,0';
cond1 = '1,2,3;1,2,3'; cond2 = '4,5,6;4,5,6'; cond3 = '7,8,9;7,8,9'; 
% go
classify_RAW_eeglab_data(filepath,filename,outpath,nFolds,cChannelSet,...
    method,crossclass,erp_baseline,cond1,cond2,cond3);

method = 'linear,dprime';


% results folder, can also be a cell array to test folder {2} against {1}
folder_name = '/Users/VU-MBP/MVPA_RESULTS_RAW/';
% settings
gsettings.mpcompcor_method = 'cluster_based' ; % can be 'none', ...
                                 % 'uncorrected', 'fdr' or 'cluster_based'
gsettings.indiv_pval = .05; % default = .05
gsettings.cluster_pval = .05; % default = .05
gsettings.one_two_tailed = 'one'; % can be 'one' or 'two'
gsettings.timelim = [-250 1500]; %  upper and lower value or empty, in ms
gsettings.timelim2 = [-250 1500]; %  upper and lower value or empty, in ms
gsettings.freqlim = [2 30]; % upper and lower value, single or empty, in Hz
gsettings.plottype = '3D'; % 2D (diagonal) or 3D (GAT or freq-time)
gsettings.channelpool = 'ALL'; % can be 'OCCIP' or other;
% go
[stats, weights, gsettings] = compute_group_MVPA(folder_name,gsettings);

% settings
gsettings.plottype = '3D'; % '2D' (line graph): accuracy * time (diagonal)
               % '3D' (color = accuracy): time * time (GAT) or freq * time
gsettings.acclim2D = [.1,.33]; % upper and lower limit in line graph
gsettings.acclim3D = [1/6-1/12 1/6+1/12]; % same for 3D hotmap color plot
gsettings.plotsigline_method = 'graphline'; 
                    % can be 'graphline' (follows the graph) or 
                    % 'straight' (below graph) or 'both' (default)
gsettings.timetick = 100; % tick used for timing
gsettings.freqtick = []; % tick used for frequency
gsettings.cent_acctick = 1/6; % chance: center accuracy (for 2D or 3D)
gsettings.acctick = .2-1/6; % tick used for accuracy (y-axis in 2D)
gsettings.downsamplefactor = 1; % downsample if plot is very spiky
gsettings.plotsubjects = false; % plot individual subjects or only average
% go
plot_MVPA(stats,gsettings);

% only one frequency, 2D
gsettings.plottype = '2D';
gsettings.freqlim = 10;

% extract values from stats
[averages, rois]= extract_MVPA_averages(stats,mask);



% topoplot
gsettings.plotweights_or_pattern = 'weight'; % can be 'covpattern' 
                                % 'corpattern' or 'weight'
gsettings.plotweights_normalized = 'normalized'; % 'non_normalized' or 
                                % 'normalized' (z-scoring across space)
gsettings.plotweights_time = 264; % time or time boundaries for which to  
                                    % plot weights in ms (takes average), 
                                    % if empty plots where accuracy peaks
gsettings.plotweights_freq = [10 14]; % indicate freq or freq boundaries  
                                      % for which to plot (takes average) 
gsettings.plotweights_weightlim = [-2.3 2.3]; % upper and lower boundary
gsettings.plotweights_imgtype = 'allpng'; % 'vec' 'png' or 'all'
% go
new_weights = plot_MVPA_weights(weights,stats,gsettings);

% tips for getting nice vector based plot:
% create a png graphic, a vec and an all
% open and image trace png in illustrator using highfidelity photo, click
% expand, ungroup, delete unnecessary bounding ring, copy to vec pdf and 
% move to back, select color bar from pdf and do object > rasterize
% (copy all to pdf to get colorized bar) and copy to vec

