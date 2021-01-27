%% set trigger labels and epoch settings
clear;
triglabels = [ 1 2 3 4 ];
files =  { 'subject01'     'subject02'    'subject03' ...
           'subject04'     'subject05'    'subject06' };

%% epoch data (e.g. through detrending)
% HERE YOU NEED SOME CODE TO EPOCH AND DETREND AND/OR HIGHPASS FILTER YOUR DATA


%% run ICA on Lisa
filepath = '/SOMEFOLDER/EEGLAB_FILES';
outpath = '/SOMEFOLDER/ICA_FILES';
compute_ICs(filepath,files,outpath,'no','no','no','yes');

%% Remove ICs using qsub
filepath = '/SOMEFOLDER/ICA_FILES';
outpath = '/SOMEFOLDER/NOBLINK_FILES';
remove_ICs(filepath,files, outpath,'0,0,1,0');
