% eeg_detrend() - Remove linear trends from epochs
%
% Usage:
%   >> EEG = eeg_detrend(EEG);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%
% Output:
%   EEG       - EEGLAB EEG structure
%
% Author: Andreas Widmann, University of Leipzig, 2006

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

% $Id: eeg_detrend.m 1 2006-04-07 12:19:45Z widmann $

function EEG = eeg_detrend(EEG)

disp('eeg_detrend(): Detrending EEG data');

EEG.data = reshape(permute(EEG.data, [2 1 3]), [EEG.pnts EEG.nbchan * EEG.trials]);
EEG.data = detrend(EEG.data);
EEG.data = permute(reshape(EEG.data, [EEG.pnts EEG.nbchan EEG.trials]), [2 1 3]);