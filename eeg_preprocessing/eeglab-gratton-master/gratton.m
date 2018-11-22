% gratton() - Applies Ocular Correction as specified by Gratton et al.,
%             1983 to the EEG-data using channel eog for regression.
%
% Usage:
%   >>  geeg = gratton( eeg, eog, blinkcritvolt, blinkcritwin );
%
% Inputs:
%   eeg     - to-be-corrected data in the format:
%             Electrode x Time x Trial
%             such that eeg(3, 400, 6) is the 3rd electrode at sampling point
%             400 of trial 6
%   eog     - data of the eog-channel -- format Time x Trial
%   blinkcritvolt - voltage sufficient for blink detection
%   blinkcritwin  - time window for blink detection (here in sampling
%                   points)
%
% Outputs:
%   ceeg    - corrected EEG data (same format as input)
%
% See also:
%   POP_GRATTON, EEGLAB

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
% $Log: gratton.m,v $
% Revision 1.1  2008/03/14 13:28:03  mihrke
% eeglab_plugins added and more stuff
%
% Revision 1.1  2007/06/22 12:29:29  mihrke
% added a lot of analysis scripts for exp 6
%
%

function ceeg = gratton( eeg, eog, crit, win );

if nargin < 2
	help gratton;
	return;
end;

[el times trials] = size(eeg);
ceeg = eeg;
assert(all([times trials] == size(eog)), 'bad input data dimensions');

avg = mean(eeg, 3); % avg trials
neeg = eeg - repmat(avg, [1 1 trials]);
neog = eog - repmat(mean(eog, 2), [1 trials]);

% detect blinks
winhalf = round(win/2);
blinks = false(size(eog));
blinks(winhalf+1:end-winhalf,:) = ((2*eog(winhalf+1:end-winhalf,:)-(eog(1:end-ceil(win),:)+eog(ceil(win)+1:end,:)))>crit);
% blinks = ((eog-circshift(eog, [winhalf 0]))+...
%     (eog-circshift(eog, [-winhalf 0])))>crit;
disp(sprintf('  gratton(): Found blinks in %i points', length(find(blinks>0))));

% loop through electrodes and get the K's
Kblinks = []; % coefficients within blinks
K = [];       % coefficients outside blinks
x = neog(:);
b = logical(blinks(:));

% correction within blinks if appropriate
if any(any(blinks))
    for e = 1:el
        y = reshape(neeg(e,:,:), [times trials]);
        y = y(:);
        Kblinks = [Kblinks regress(y(b), x(b))];
    end;
    disp(' gratton(): Coefficients within blinks:');
    disp(Kblinks);
    minus = (repmat(Kblinks', [1 times trials]).*shiftdim(repmat(eog, [1 1 el]), 2));
    ceeg(blinks) = eeg(blinks) - minus(blinks);
    avg = mean(ceeg, 3); % avg trials
    neeg = ceeg - repmat(avg, [1 1 trials]);
    disp(' gratton(): Corrected EEG within blinks');
end;

% correction outside of blinks
for e = 1:el
    y = reshape(neeg(e,:,:), [times trials]);
    y = y(:);
    K = [K regress(y, x)];
end;

disp(' gratton(): Coefficients outside of blinks:');
disp(K);
ceeg = ceeg - (repmat(K', [1 times trials]).*shiftdim(repmat(eog, [1 1 el]), 2));
disp(' gratton(): Corrected EEG outside of blinks');

% figure;
% for i = 1:times
%     plot(1:times, eog(:,i));
%     hold on;
%     %plot(1:times, neeg(1, :, 1), 'r');
%     plot(1:times, blinks(:,i)*50, 'k');
%     hold off;
%     x = waitforbuttonpress;
% end;

return;
