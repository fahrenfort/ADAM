function output = file_list_grouptrain(listA,listB)
% file_list_grouptrain(listA,listB)
% listA and listB are cell arrays of strings, in which listA contains the train files for all
% subjects and listB contains the test files for all subjects
% the output list contains semicolon separated combinations of all train and test files, but without
% ever using the same subject for training and testing
% This function works under the assumption that train and test files of the same subject are
% identifiable using a unique number in the filename (e.g. 'subj01').
% 
%
% J.J.Fahrenfort, UvA/VU 2018

cOut = 0;
for cA = 1:numel(listA)
    for cB = 1:numel(listB)
        %subjA = regexp(listA{cA},'*\d*','Match');
        %subjB = regexp(listB{cB},'*\d*','Match');
        subjA = listA{cA};
        subjB = listB{cB};
        % if neither of the strings is contained in the other, or if the first element of the
        % testfile is a digit when removing the part of the pattern that is present in the trainfile
        % (e.g. when trainfile is 'subj1' and testfile is 'subj12')
        if (isempty(strfind(subjA,subjB)) && isempty(strfind(subjB,subjA))) || ~isempty(find(isstrprop(regexprep(subjB,subjA,''),'digit'),true,'first')==1)
            cOut = cOut + 1;
            output{cOut} = [listA{cA} ';' listB{cB}];
        end
    end
end
output = sort(output);