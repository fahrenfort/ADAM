function varargout = funname(varargin)
% This function is a backward compatibility wrapper based on Robert Oostenveld's
% wrappers for FieldTrip.
%
% Please type 'help nt_xxx' where xxx is the name of you just typed.
% 
% NoiseTools.

eval(['persistent ',mfilename,'_been_here_before;']);

if eval(['isempty(',mfilename,'_been_here_before)'])
disp(['Warning: ', mfilename, ' is an old-style name, replace by nt_', mfilename]);
eval([mfilename,'_been_here_before=1;']);
end

prefix    = 'nt_';
funname   = mfilename;
funhandle = str2func([prefix funname]);
[varargout{1:nargout}] = funhandle(varargin{:});
