function output = file_list_restrict(varargin)
% function output = file_list_restrict(varargin)
% filter a cell array of files by one or more criteria
% first argument is a cell array of file names
% all ensuing arguments can be strings that are used to select from the
% file names.
% Examples:
% condA_files = file_list({'subj1_condA', 'subj1_condB', subj2_condA', subj2_condB'},'condA');
% subj2_condB_files = file_list({'subj1_condA', 'subj1_condB', subj2_condA', subj2_condB'},'condB','subj2');
%
% By J.J.Fahrenfort, VU, 2015

output = varargin{1};
for c = 2:numel(varargin)
    ind2remove = [];
    for cFiles = 1:numel(output)
        if isempty(strfind(output{cFiles},varargin{c}))
            ind2remove = [ ind2remove cFiles ];
        end
    end
    output(ind2remove) = [];
end