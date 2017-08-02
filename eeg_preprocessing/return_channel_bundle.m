function [channelset, bundlename_or_bundlelabels] = return_channel_bundle(channelset)
% internal ADAM function to determine the channel bundle name and/or cell
% array of the electrodes that were specified in channelset
%
% input can be:
% channelset  -  a bundle name as a string, a csv string containing electrodes, or a bundle number 
%
% output is:
% channelset  -  either a bundle name as string or a string of concatenated electrodes as a string separated by underscores 
% bundlename_or_bundlelabels  -  either a bundle name as string or a cell array of electrodes
if isempty(channelset)
    channelset = 1;
end
if ischar(channelset)
    bundlename_or_bundlelabels = regexp(channelset,',','split'); % now a cell array of strings
    if numel(bundlename_or_bundlelabels) == 1
        bundlename_or_bundlelabels = bundlename_or_bundlelabels{1}; % now a string
    end
    if iscell(bundlename_or_bundlelabels) % a cell array of electrodes
        channelset = regexprep(channelset,',','_'); % a string
    elseif ~isnan(string2double(channelset))
        channelset = string2double(channelset); % a number
    end
end
if isnumeric(channelset) % for backward compatibility, in principle channelset should be a single string or csv of electrodes
    if channelset > 0
        bundlenames = {'ALL' 'OCCIP' 'PARIET' 'FRONTAL' 'TEMPORAL' 'OCCIPARIET' 'CDA' 'N2Pc_SPCN' };
        channelset = bundlenames{channelset}; % now a string
    else
        channelset = 'ALL_NOSELECTION'; % now a string
    end
    bundlename_or_bundlelabels = channelset; % is now also a string
end