function [condname, suborder] = find_plot_order(dirz,cdirz,graphsettings)
% find where to put the plot and what name to put above the plot

thisdir = dirz{cdirz};
condname = strrep(thisdir,'_',' ');
% find the plot order from the order of the text patterns in graphsettings.plotorder
if exist('graphsettings') && isfield(graphsettings,'plotorder') && ~isempty(graphsettings.plotorder)
    suborder = find(strcmp(repmat({thisdir},1,numel(graphsettings.plotorder)),graphsettings.plotorder));
    if isempty(suborder)
        suborder = find(~cellfun(@isempty,(cellfun(@strfind,repmat({thisdir},1,numel(graphsettings.plotorder)),graphsettings.plotorder,'UniformOutput',false))),1);
    end
    if isempty(suborder)
        error(['cannot find a pattern in graphsettings.plotorder that matches ' thisdir '. Specify graphsettings.plotorder = [] or use patterns that match your conditions.' ]); 
    end
else
    suborder = cdirz;
end