function common_to_use = commonstring(Scell)
% Find the part of the start of a string that is common in a cell array of
% strings. Returns empty if there is no common part.
common_to_use = Scell{1}(1:find(any(diff(char(Scell(:))),1),1,'first')-1);
if isempty(common_to_use)
   common_to_use ='';
end