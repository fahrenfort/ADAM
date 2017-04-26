function [x, reg] = gsmooth(x, sigma, varargin)
%GSMOOTH Gaussian vector smoothing
%   [SMOOTH, REGION] = GSMOOTH(X, SIGMA) carries out Gaussian smoothing on
%   a sequence.
%
%   X must be a vector of class double or single.
%
%   SIGMA specifies the smoothing constant (i.e. the "sigma" parameter of
%   the Gaussian smoothing kernel).
%
%       The output vector will be smaller than the input vector in order to
%       avoid having to extrapolate beyond the ends of the vector. The
%       result REGION reports the range of indices in the input for which
%       corresponding output values have been computed, in the form [MINI,
%       MAXI]. The length of the output vector is MAXI-MINI+1. The
%       reduction in size depends on the smoothing constant, and is chosen
%       so that the smoothing mask is a good approximation to the
%       untruncated Gaussian mask.
%
%   REGION = GSMOOTH(X, SIGMA, 'RegionOnly') returns only the region.
%
%   [SMOOTH, REGION] = GSMOOTH(X, SIGMA, NAME, VALUE, ...) allows
%   additional parameters to be set using name-value pairs. Names may be as
%   follows:
%
%       'Region' - Allows a region of interest in X to be specified. Its
%       value may be:
%
%           'valid' or [] - (default). As described above.
%
%           'same' - The region of interest covers the whole of X, and the
%           output vector has the same size as X. Reflection at the
%           boundaries is used to extrapolate X if it does not wrap around.
%
%           A 2-element row vector with elements [MINI, MAXI] giving the
%           minimum and maximum indices of a segment of X. The results
%           vector will have length MAXI-MINI+1 and will contain the
%           smoothed values for the specified region of the input vector.
%           Reflection at the ends is used for extrapolation of X when
%           necessary.
%
%       'Wrap' - Specifies that the vector wraps round (i.e. has periodic
%       boundary conditions). The value must be a logical scalar. The
%       default is false. If the vector wraps, the default output region
%       will be the same size as the vector. [] may be given as an
%       alternative to false.
%
% See also: GSMOOTH2, GSMOOTHN

% Copyright David Young 2014

if ~isvector(x)
    error('DavidYoung:gsmooth:notVector', ...
        'First argument must be a vector');
end

row = isrow(x);
if row
    x = x.';
end

if ~isempty(varargin) && strcmp(varargin{1}, 'RegionOnly')
    % deal with region only case
    x = gsmoothn(x, sigma, varargin{:});
else
    [x, reg] = gsmoothn(x, sigma, varargin{:});
    if row
        x = x.';
    end
end

end

