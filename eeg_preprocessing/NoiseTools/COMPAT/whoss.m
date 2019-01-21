function varargout = funname(varargin)

error('use nt_whoss instead');
return

% code below doesn't work for whoss because nt_whoss doesn't have access to
% caller workspace

% This function is a backward compatibility wrapper based on Robert Oostenveld's
% wrappers for FieldTrip.
%
% Please look in nt_xxx for the help of the function that you are looking
% for, where xxx is the name of the function that you were looking for.
eval(['persistent ',mfilename,'_been_here_before;']);

if eval(['isempty(',mfilename,'_been_here_before)'])
disp(['Warning: ', mfilename, ' is an old-style name, replace by nt_', mfilename]);
eval([mfilename,'_been_here_before=1;']);
end

prefix    = 'nt_';
funname   = mfilename;
funhandle = str2func([prefix funname]);
[varargout{1:nargout}] = funhandle(varargin{:});
