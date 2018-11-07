function remove_ICs(filepath,filename, outpath, hor_ver_blink_gd)
% function remove_ICs(filepath,filename, outpath, hor_ver_blink_gd)
% removes ICAs from epoched data.
% Uses the ADJUST plugin and/or a text file to specify which components to
% remove. hor_ver_blink_gd determines whether you want to remove hor
% (horizontal eye movements), ver (vertical eye movements), blinks, and/or
% gd (generic discontinuities). For example, specify hor_ver_blink_gd as:
% '1,1,1,0' to remove everything except generic discontinuities.
% Removes all artifactual ICs by default if left empty or if it is not
% correctly specified.
% It is possible to input a different sourcefile or sourcedirectory for
% determining which components to remove. This is done by appending a comma
% and a second filename or foldername after filename or filepath. The
% second file/folder name will be the source file of the ICA. E.g.: use
% 'PP01test,PP01train' to use the ICA from PP01train, and remove the
% identified components in PP01test. The ICA matrix from the source file
% will be copied over and applied to the destination file.
% Finally, it is possible to indicate components to remove by
% modifying a text file in the ICA or destination folder (ICA folder takes
% precedence). This text file has the exact same name as the filename, but
% appended with .txt, containing components to remove separated by single
% spaces, e.g.: 2 3 45 1. Initially, this file is created and filled with
% the components identified by ADJUST, which you can then modify after
% visual inspection and re-run the function, after which it will only
% remove the indicated components. If left empty, no components will be
% removed.
%
% Example usage:
% remove_ICs_and_epoch(pwd,'ICA_WM_pilot_pp1',pwd,hor_ver_blink_gd);
% create_qsub_files('$HOME/Dropbox/matlab_scripts/eeg_preprocessing','remove_ICs',settings,'$HOME/EEG_LAB_DATA/ICA',{'WM_PP1,WM_train_PP1','WM_PP2,WM_train_PP2'},'$HOME/EEG_LAB_DATA/CLEAN'));
%
% By J.J.Fahrenfort, VU 2014, 2016

if isempty(outpath)
    outpath = filepath;
end
if ~exist(outpath,'dir')
    mkdir(outpath);
end
if nargin < 4
    hor_ver_blink_gd = [];
end
if ischar(hor_ver_blink_gd)
    hor_ver_blink_gd = str2num(hor_ver_blink_gd);
end
if isempty(hor_ver_blink_gd) || numel(hor_ver_blink_gd)~=4
    disp('removing only EOG ICs (default)');
    hor_ver_blink_gd = [1 1 1 0];
end

% get filepaths
filepaths = regexp(filepath,',','split');
filepath = filepaths{1};
if numel(filepaths) > 1
    icapath = filepaths{2};
else
    icapath = filepath;
end
if ~exist(filepath,'dir')
    error('the input directory does not exist');
end

% get filenames
filenames = regexp(filename,',','split');
[~,filename1,~] = fileparts(filenames{1});
if numel(filenames) > 1
    [~,filename2,~] = fileparts(filenames{2});
else
    filename2 = filename1;
end

% ICs2reject
ICs2rejectFile = [icapath filesep filename2 '.txt'];
ICsrejectedFile = [outpath filesep filename1 '.txt'];
if ~exist(ICs2rejectFile,'file') && exist(ICsrejectedFile,'file')
    ICs2rejectFile = ICsrejectedFile;
end
% bookkeeping
delete([outpath filesep filename1 '*comp*.png']);

fprintf(['\n\nRemoving components from EEG file ' filepath filesep filename1 ]);
fprintf(['\nUsing ICA from source file ' icapath filesep filename2 '\n\n']);

% load destination EEG (the file on which to work)
EEG = pop_loadset('filename',[filename1 '.set'],'filepath',filepath);

% make sure it contains channel information
EEG = add_channel_info(EEG);

% load source EEG if present (the file from which to take ICA)
EEG2 = pop_loadset('filename',[filename2 '.set'],'filepath',icapath);
EEG2 = add_channel_info(EEG2);

% making sure both sets have the same channels
channelnames = {EEG.chanlocs(:).labels};
channelnames2 = {EEG2.chanlocs(EEG2.icachansind).labels}; % only selecting channels for which ICA was obtained
channels2keep = intersect(channelnames,channelnames2);
EEG = pop_select(EEG,'channel',channels2keep);
EEG2 = pop_select(EEG2,'channel',channels2keep);
if sum(strcmpi({EEG.chanlocs(:).labels},{EEG2.chanlocs(:).labels})) ~= EEG.nbchan
    error('The electrode order should be the same in both datasets! Needs fixing.');
end
disp(['using the following channels:  ' sprintf('%s, ',channels2keep{1:end-1}) channels2keep{end}]);

% copy over source
disp('copying over ICA info from source');
EEG.icawinv = EEG2.icawinv;
EEG.icasphere = EEG2.icasphere;
EEG.icaweights = EEG2.icaweights;
EEG.icachansind = EEG2.icachansind;

% recompute icaact
disp('Recomputing EEG.icaact from EEG.icaweights, EEG.icasphere and EEG.data');
EEG.icaact = reshape(EEG.icaweights*EEG.icasphere*reshape(EEG.data(1:size(EEG.icaweights,1),:,:),[size(EEG.icaweights,1) size(EEG.data,2)*size(EEG.data,3)]),[size(EEG.icaweights,1) size(EEG.data,2) size(EEG.data,3)]);

% either retrieve components from file or use automatic rejection
if exist(ICs2rejectFile,'file')
    disp(['reading components to remove from text file ' ICs2rejectFile]);
    ICs2reject = textread(ICs2rejectFile);
else
    % run adjust on source, only works when destination data is epoched!
    % (ICA source data can be continuous if desired)
    [art, horiz, vert, blink, disc ] = ADJUST(EEG,[outpath filesep 'ADJUST_' filename1 '.txt' ]);
    disp('using components that were just computed');
    ICs2reject = unique([horiz*hor_ver_blink_gd(1) vert*hor_ver_blink_gd(2) blink*hor_ver_blink_gd(3) disc*hor_ver_blink_gd(4)]);
end
ICs2reject = ICs2reject(ICs2reject>0);

% remove ICA components from destination and write removed components to file
if isempty(ICs2reject)
    disp('no components were identified/removed');
else
    % save property plots of to be removed ICs (do not include EOG, for clarity)
    for c=1:numel(ICs2reject)
        disp(['removing component: ' num2str(ICs2reject(c))]);
        pop_prop_savepng(pop_select(EEG2, 'channel', select_channels({EEG2.chanlocs(:).labels},'EEG')), 0, ICs2reject(c), NaN, {'freqrange' [2 50] },[outpath filesep filename2 '_srce']);
        pop_prop_savepng(pop_select(EEG, 'channel', select_channels({EEG.chanlocs(:).labels},'EEG')), 0, ICs2reject(c), NaN, {'freqrange' [2 50] },[outpath filesep filename1 '_dest']);
    end
    clear EEG2;
    % and remove
    EEG = pop_subcomp(EEG, ICs2reject, 0);
end
fid = fopen(ICsrejectedFile,'w');
fprintf(fid,'%d ',ICs2reject);
fclose(fid);

% write results
pop_saveset(EEG, 'filename',[filename1 '.set'],'filepath',outpath);

