% Copyright (C) 2007  Matthias Mittner <matthias.mittner@uit.no>
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% $Log: eegplugin_gratton.m,v $
% Revision 1.1  2008/03/14 13:28:03  mihrke
% eeglab_plugins added and more stuff
%
% Revision 1.1  2007/06/22 12:29:29  mihrke
% added a lot of analysis scripts for exp 6
%
%
function eegplugin_gratton( fig, try_strings, catch_strings);

% create menu
toolsmenu = findobj(fig, 'tag', 'tools');

% build command for menu callback
cmd = [ 'EEG = pop_gratton(EEG); [ALLEEG EEG CURRENTSET]=eeg_store(ALLEEG,EEG);'];
cmd = [ cmd 'eeglab redraw;' ];

finalcmd = [ try_strings.no_check cmd ];
finalcmd = [ finalcmd 'LASTCOM = ''' cmd ''';' ];
finalcmd = [ finalcmd catch_strings.store_and_hist ];

submenu = uimenu( toolsmenu, 'label', 'Ocular Correction (Gratton)',...
    'callback', finalcmd);
