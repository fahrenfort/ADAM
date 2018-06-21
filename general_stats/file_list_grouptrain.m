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
        subjA = regexp(listA{cA},'\d*','Match');
        subjB = regexp(listB{cB},'\d*','Match');
        if ~strcmp(subjA,subjB)
            cOut = cOut + 1;
            output{cOut} = [listA{cA} ';' listB{cB}];
        end
    end
end