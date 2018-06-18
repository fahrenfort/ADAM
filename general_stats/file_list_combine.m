function output = file_list_combine(listA,listB,patternlist)
% function output =  file_list_combine(listA,listB,patternlist)
% listA and listB are cell arrays of strings
% patternlist is a cell array of criteria that are used to combine listA and
% listB, where each new cell separates the matched files from listA and
% listB by a comma and puts them in a new cell. If patternlist is omitted,
% it simply matches the elements one by one from the ordered lists.
% Example:
% new_list = file_list_combine({'subj1_condA', 'subj2_condA'},{'subj1_condB', 'subj2_condB'},{'subj1','subj2'});
% or
% new_list = file_list_combine({'subj1_condA', 'subj2_condA'},{'subj1_condB', 'subj2_condB'});
%
% By J.J.Fahrenfort, VU, 2015, 2016
if nargin > 2
    output{numel(patternlist)} = [];
    for c = 1:numel(patternlist)
        indA = [];
        indB = [];
        for cA = 1:numel(listA)
            if ~isempty(strfind(listA{cA},patternlist{c}))
                indA = cA;
            end
        end
        for cB = 1:numel(listB)
            if ~isempty(strfind(listB{cB},patternlist{c}))
                indB = cB;
            end
        end
        output{c} = [listA{indA} ';' listB{indB}];
    end
else
    if numel(listA) ~= numel(listB)
        error('unequal number of elements in list A and list B');
    else
        listA = sort(listA);
        listB = sort(listB);
    end
    output{numel(listA)} = [];
    for c = 1:numel(listA)
        output{c} = [listA{c} ';' listB{c}];
    end
end