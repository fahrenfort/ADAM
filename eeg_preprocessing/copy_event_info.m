function copy_event_info(filepaths,filenames)
% function copy_ICA_info(filepaths,filenames)
% copies the event information in one set of EEGlab files to a different
% set of EEGlab files.
% filenames is assumed to be a cell array of file names.
% If the file names are different for files, separate them by a
% comma (first = source, second = target), otherwise the same filename
% names are assumed for source and target.
% If the filepaths are different for both files, separate them by a
% comma (first = source, second = target), otherwise it is assumed all
% files are in the same directory.
% Example usage:
% copy_event_info('/Users/myname/data/infiles,/Users/myname/data/outfiles',{'subj1_file'},{'subj_file'});
% copy_event_info([pwd filesep 'source,' pwd filesep 'target'],{'subj1_file1,subj1_file2'},{'subj2_file1,subj2_file2'});
% Johannes Fahrenfort, VU 2016

filepaths = regexp(filepaths,',','split');
inpath = filepaths{1};
if numel(filepaths)>1
    outpath = filepaths{2};
else
    outpath = inpath;
end

if nargin<2
    % no filenames specified, taking all files in the input directory
    filenames = dir(fullfile(inpath,'*.set'));
    filenames = {filenames(:).name};
end

for files = filenames
    fnames = regexp(files{1},',','split');
    [~,infile] = fileparts(fnames{1});
    infile = [infile '.set'];
    if numel(fnames) > 1
        [~,outfile] = fileparts(fnames{2});
        outfile = [outfile '.set'];
    else
        outfile = infile;
    end
    if ~exist(fullfile(inpath,infile),'file')
        error([fullfile(inpath,infile) ' does not exist.']);
    end
    if ~exist(fullfile(outpath,outfile),'file')
        error([fullfile(outpath,outfile) ' does not exist.']);
    end
    if strcmp(fullfile(inpath,infile),fullfile(outpath,outfile))
        error(['source and target are the same:' fullfile(inpath,infile)]);
    end
    disp(['Loading subject ' infile ]);
    EEG = pop_loadset('filename',infile,'filepath',inpath);
    % getting event values
    eventstruct = EEG.event;
    if isfield(EEG,'epoch')
        epochstruct = EEG.epoch;
    end
    EEG = pop_loadset('filename',outfile,'filepath',outpath);
    % copying over to new file
    EEG.event = eventstruct;
    EEG.epoch = epochstruct;
    EEG = eeg_checkset(EEG);
    % saving result
    pop_saveset(EEG, 'filename',outfile,'filepath',outpath);
end
disp('done!');
