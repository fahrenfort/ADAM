function [ correct, outcomes ] = checkvalues(bdfvalue,newvalue)
% function to check whether new triggervalues make sense given wrong
% bitcodes in the bdf file
binNEW = dec2bin(newvalue,8);
new2bdf = bin2dec(binNEW(:,end-6:end));
if sum(size(new2bdf) ~= size(bdfvalue))
    new2bdf = new2bdf';
end
outcomes = bdfvalue == new2bdf;
correct = sum(outcomes) == numel(bdfvalue);