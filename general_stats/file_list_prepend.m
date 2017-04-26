function output = file_list_prepend(file_list,string)
% function output = file_list_prepend(file_list,string)
% Prepend a string to each cell in a cell array
% By J.J.Fahrenfort, VU, 2015

output{numel(file_list)} = [];
for c = 1:numel(file_list)
    output{c} = [string file_list{c}];
end