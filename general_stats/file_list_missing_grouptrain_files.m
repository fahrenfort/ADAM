function missing_files = file_list_missing_grouptrain_files(grouptrainfiles,folder_name)
% Determines whether all output files were generated in the results folder, and if not which files
% are missing. When called without specifying a resultsfolder the opens a selection dialog asking to
% select the folder to which this pertains. As output generates a list of the missing train-test
% combinations.
% 
% Inputs:
%       grouptrainfiles     =   a cell array containing the train test file combinations that were
%                               entered into the analysis. E.g. cfg.grouptrainfiles =
%                               {'subj1train,subj2test', 'subj2train, subj1test'}; The function
%                               will determine whether each grouptrain file generated a results
%                               file
%       folder_name         =   the folder for which the missing grouptrain files are checked.
%
if nargin < 2
    folder_name = '';
end

if isempty(folder_name)
    folder_name = uigetdir('','select results directory for which to determine missing files');
    if ~ischar(folder_name)
        error('no folder was selected');
    end
end

allfiles = dir([folder_name filesep '*.mat']);
allfiles = {allfiles(:).name};

for c=1:numel(allfiles)
    file = allfiles{c};
    trainfile = file(strfind(file,'train')+6:strfind(file,'test')-2);
    testfile = file(strfind(file,'test')+5:strfind(file,'.')-1);
    allfiles{c} = ([trainfile ';' testfile]);
end
missing_files = grouptrainfiles(~ismember(grouptrainfiles,allfiles));
disp(missing_files');
