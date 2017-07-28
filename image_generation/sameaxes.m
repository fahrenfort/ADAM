function sameaxes(XYZC, HNDS)
%SAMEAXES unifies/synchronizes axis limits on different axes and subplots.
% Calling sameaxes() sets identical limits [the pooled min() and max()] to
% each axis respectively - which is very handy to compare different plots!
% It is possible to restrict the effect to a specific figure and/or axis:
% Use XYZC to select which axis to unify (e.g. only the y- or color-axis),
% and HNDS to restrict the operations to  e.g. only the current figure.
% Note that it doesn't link any axes as done by linkaxes().
% You may exclude a figure or axis by setting 'HandleVisibility' to 'off'.
% Many examples below.
%
% INPUT
%   XYZC - which axis to sync (cell or char array) [default 'xyzc' ==> all]
%   HNDS - figure handle(s) to search for children [default 0 ==> all axes]
%
% OUTPUT
%   none; re-sets the [xyzc]lim-properties of HNDS objects' children axes
%
% EXAMPLES
%   sameaxes() without arguments unifies ALL existing axes (x, y, z, color)
%   sameaxes('y') applies to the ylim of ALL existing axes (across figures)
%   sameaxes('xc', gcf()) unifies only the x-axis- and color-limits of all
%                         subplots that are children of the current figure
%   sameaxes([], [fg1,fg2]) separately unifies children of figures fg1, fg2
%
%   In case you want to exclude certain figures or axes from being
%   affected, you can set their handles' visibilities to 'off' beforehand.
%   Here an example with 5 figures, excluding the first and the fifth:
%   for ii = 1:5
%       figure(), plot(ii*rand(10)), title(sprintf('plot %d', ii))
%   end
%   figures = flipud(findobj('Type','figure')); % get figure handles
%   exclude = figures([1 5]);                   % select 1,5 for exclusion
%   set(exclude, 'HandleVisibility', 'off')     % set handles "invisible"
%   sameaxes()                                  % sync all visible handles
%   set(exclude, 'HandleVisibility', 'on')      % restore to default
%
% See also PBASPECT, DASPECT, LINKAXES, XLIM, YLIM, ZLIM.

% Created Jun/13 by Johannes Keyser (jkeyser@uni-osnabrueck.de)
% Revised Aug/13 by jkeyser: +arg HNDS to restrict axes search from parents
% Revised Okt/13 by jkeyser: +arg XYZC to generalize to any axis, +comments
% Revised Jan/14 by jkeyser: +polished for Matlab File Exchange publication
% Revised May/14 by jkeyser: +fixed check for handles, +exclusion example

validXYZC = 'xyzc'; % x-, y-, z- and/or color-axes (in any combination)
if nargin < 1 || isempty(XYZC), XYZC = validXYZC; end
if nargin < 2 || isempty(HNDS), HNDS = 0;         end % 0 is root object
% parse XYZC input regardless of cell vs char (e.g. {'x','y','c'} vs 'xyz')
if iscell(XYZC)
    XYZC = cellfun(@(c) lower(c), XYZC, 'UniformOutput', false);
    XYZC = char([XYZC{:}]);
elseif ischar(XYZC)
    XYZC = lower(XYZC);
end
assert(ischar(XYZC) && all(arrayfun(@(c) ismember(c,validXYZC), XYZC)),...
    'XYZC must be a subset, or combination of {''x'',''y'',''z'',''c''}!')
assert(all(ishandle(HNDS)), 'Argument HNDS must be valid object handles!')
for hnd = HNDS(:)' % iterate over parent handles [sync parents separately!]
    for xyz = XYZC(:)' % iterate over x, y, z, c
        lim = [xyz 'lim'];
        % find axes (== objects with color-limits - except colorbars)
        axs = findobj(HNDS, '-property',lim, '-not','tag','colorbar');
        if isempty(axs)
            warning('sameaxes:noSuchLim', 'No children with "%s".', lim)
            continue
        end
        % get() and then set() the pooled min() & max() to all axes objects
        lims = get(axs, lim);
        if iscell(lims)
            lims = [lims{:}];
        end % unpack if necessary
        set(axs, lim, [min(lims) max(lims)])
    end
end
end