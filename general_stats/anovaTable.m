function [s] = anovaTable(AT)
c = table2cell(AT);
% The code below replaces the blank(?) entries in the F and p columns
% with empty character strings.  Not sure why, but the blank entries are
% actually 1 (F) or 0.05 (p).  
%
% The code below can probably be simplified, but I'm not sure exactly how.
%
% An assumption here is that measurements are made on human subjects
% ("Participant").  If the measurements are on some other entity, such as
% "flowers" in the example in the ranova documentation, then Participant,
% below, should be changed accordingly.
for i=1:size(c,1)       
        if c{i,4} == 1
            c(i,4) = {''};
        end
        if c{i,5} == .5
            c(i,5) = {''};
        end
end
effect = AT.Properties.RowNames;
for i=1:length(effect)
    tmp = effect{i};
    tmp = erase(tmp, '(Intercept):');
    tmp = strrep(tmp, 'Error', 'Participant');
    effect(i) = {tmp}; 
end
% determine the required width of the Effects column
fieldWidth = max(cellfun('length', effect));
barDouble = sprintf('%s\n', repmat('=', 1, fieldWidth + 57));
barSingle = sprintf('%s\n', repmat('-', 1, fieldWidth + 57));
% re-organize the data 
c = c(2:end,[2 1 3 4 5]);
c = [num2cell(repmat(fieldWidth, size(c,1), 1)), effect(2:end), c]';
% create the ANOVA table
s = barDouble;
s = [s sprintf('%-*s %4s %10s %14s %10s %10s\n', fieldWidth, 'Effect', 'df', 'SS', 'MS', 'F', 'p')];
s = [s barSingle];
%s = [s, sprintf('%-*s %4d %14.3f %14.3f %10.3f %10.4f\n', c{:})];
s = [s, sprintf('%-*s %4d %14.3f %14.3f %10.3f %.2e\n', c{:})];
s = [s, barDouble];
end