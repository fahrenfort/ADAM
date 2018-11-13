function varargout = funname(varargin)
% This function is a backward compatibility wrapper based on Robert Oostenveld's
% wrappers for FieldTrip.
%
% Please look in nt_xxx for the help of the function that you are looking
% for, where xxx is the name of the function that you were looking for.

eval(['persistent ',mfilename,'_been_here_before;']);

if eval(['isempty(',mfilename,'_been_here_before)'])
disp(['Warning: nt_colorlines is an old-style name, replace by nt_line_colors']);
eval([mfilename,'_been_here_before=1;']);
end

[varargout{1:nargout}] = nt_linecolors(varargin{:});
