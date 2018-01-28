function FT_EEG = fix_dimord(FT_EEG,new_dimord,field)
% fix_dimord fixes the order of dimensions in FT_EEG to channel x time x trial (default) or another
% order if desired. 
% ADAM toolbox internal function by J.J.Fahrenfort
if nargin<3
    field = 'trial';
end
if nargin<2
    new_dimord = 'chan_time_rpt';
end

% if permute
if ~strcmp(FT_EEG.dimord,new_dimord)
    % original dimensions
    dims = regexp(FT_EEG.dimord, '_', 'split');
    % new dimensions
    new_dims = regexp(new_dimord, '_', 'split');
    % find order
    new_ind_dimord = [find(strcmp(dims,new_dims{1})) find(strcmp(dims,new_dims{2})) find(strcmp(dims,new_dims{3}))];
    % permute
    FT_EEG.(field) = permute(FT_EEG.(field),new_ind_dimord);
    FT_EEG.dimord = new_dimord;
end