function merge_EEG_sessions(filepath,filenames,outpath,subjects,reject_artifacts)
% function merge_EEG_sessions(filepath,filenames,outpath,subjects,reject_artifacts)
% concatenates all sessions from EEG lab file belonging together
% filepath: where files are stored
% filenames: names of files either in a cell array or as comma separated
% string. Wildcards * and ? can be used, e.g. filenames = '*.set' will
% take all the .set files in the input filepath as sources for merging
% subjects. When  using together with qsub, use '*,*.set' to specify the
% filenames to prevent qsub from expanding to a separate line  for each
% filename. Set settings.use_scratch = false; in this case. Files are
% appended in alphabetical order when using a wildcard. When a specific
% order is required for merging, specify each file separately without
% wildcards (in the correct order of course).
% Subjects: specify regular expressions for which the
% filenames will be merged. This patterns is also used for the new
% filename, after removal of the asterisk (*), dot (.) and question mark
% (?) characters, e.g. '.*subj02.*_AB.*' as a subject name will result in
% the selection of all files with that data pattern, and the file name in
% this case would be 'subj02_AB'. You can find more about how regular
% expressions work by typing: help regexp in the command line.
% NOTE: REGULAR EXPRESSIONS DO NOT USE * AS A WILDCARD, RATHER SAY .* TO
% INDICATE A WILDCARD (. means any character and * means 0 or more times)
% If you simply do 'subj02' without anything else, it will grab all files
% containing subj02 anywhere and merge them, no wildcards necessary.
% If reject_artifacts = true (default
% false), the artifacts that are contained in EEG.reject.rejmanual are
% rejected before concatenating.
% example:
% merge_EEG_sessions('c:\inputfiles','*.set','c:\mergedfiles',{'.*subj02.*_AB.*' '.*subj03.*_AB.*' '.*subj04.*_AB.*'},true);
%
% J.J.Fahrenfort, VU, 2015

% some input checking
if nargin<5
    reject_artifacts = false;
end
if nargin<4
    outpath = pwd;
end
if ischar(reject_artifacts)
    if strcmp(reject_artifacts,'true') || str2double(reject_artifacts) == 1
        reject_artifacts = true;
    else
        reject_artifacts = false;
    end
end
if ~exist(outpath,'dir')
    mkdir(outpath);
end
if ~isempty(strfind(filenames,'*')) || ~isempty(strfind(filenames,'?'))
    if ~strcmp(filenames(end-3:end),'.set')
        filenames = [filenames '.set'];
    end
    if ~isempty(strfind(filenames,','))
        filenames = filenames(max(strfind(filenames,','))+1:end);
    end
    filenames = dir([filepath filesep filenames]);
    filenames = filenames(~startsWith({filenames.name}, '.'));
    filenames = {filenames(:).name};
    % make sure file names are sorted in alphabetical order in this case
    filenames = sort(filenames);
end
if ~iscell(filenames)
    filenames = regexp(filenames, ',', 'split');
end
if ~iscell(subjects)
    subjects = regexp(subjects, ',', 'split');
end
% go
for cSubj = 1:numel(subjects)
    ALLEEG = []; EEG=[]; CURRENTSET=[];
    subjIndex = [];
    for c = 1:numel(filenames)
        if ~isempty(regexp(filenames{c},subjects{cSubj}))
            subjIndex = [subjIndex c];
        end
    end
    for cIndex = 1:numel(subjIndex)
        [~,filename,~] = fileparts(filenames{subjIndex(cIndex)});
        % here the specifics start
        EEG = pop_loadset('filename',[filename '.set'],'filepath',filepath);
        disp([filename '.set']);
        % remove trials marked by variable rejected_trials
        if reject_artifacts
            EEG = pop_select(EEG,'notrial',find(EEG.reject.rejmanual));
        end
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    end
    OUTEEG = pop_mergeset( ALLEEG, 1:numel(subjIndex), 1);
    filename = subjects{cSubj};
    filename(regexp(filename,'[*.?[]^()]')) = [];
    filename = [filename '_merged'];
    disp(['the above files were merged into ' filename '.set']);
    pop_saveset(OUTEEG, 'filename',[filename '.set'],'filepath',outpath);    
end
