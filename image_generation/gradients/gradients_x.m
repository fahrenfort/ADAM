function [xg, region] = gradients_x(x, varargin)
%GRADIENTS_X Estimate gradient of a sequence with Gaussian smoothing
%   [VG, REGION] = GRADIENTS_X(V) estimates the gradient of V using centred
%   differences.
% 
%       V is a vector of class double or single.
%
%       VG is a vector of estimates of the gradient, computed using
%       symmetric differencing between the two nearest neighbours of each
%       element of V; VG(I) = (V(I+1)-V(I-1))/2.
% 
%       REGION is the region of V for which gradients have been computed,
%       given as a 2-vector with the minimum and maximum indices. In the
%       non-smoothing case it will be [2 length(V)-1].
% 
%   [VG, REGION] = GRADIENTS_X(V, SIGMA) carries out Gaussian smoothing and
%   differencing to estimate the gradient of V.
%
%       SIGMA specifies the smoothing constant. It must be a scalar giving
%       the "sigma" parameter of a Gaussian mask. If SIGMA is zero, no
%       smoothing is done.
% 
%       The output vector will be smaller than the input vector in order to
%       avoid having to extrapolate beyond its ends. The result REGION
%       reports the region of the input vector for which the output values
%       have been estimated, in the form [MINI, MAXI] where MINI and MAXI
%       are the first and last indices of the region. The length of the
%       output vector is MAXI-MINI+1. The reduction in size depends on the
%       smoothing parameter, and is chosen so that the smoothing mask is a
%       good approximation to the untruncated Gaussian mask.
%
%   REGION = GRADIENTS_X(V, SIGMA, 'RegionOnly') returns only the region.
% 
%   [VG, REGION] = GRADIENTS_X(V, SIGMA, NAME, VALUE, ...)
%   [VG, REGION] = GRADIENTS_X(V, NAME, VALUE, ...) 
%   set additional parameters using name-value pairs. Names may be:
% 
%       'Region' - Specifies a region of interest in V. Its value may be:
% 
%           'valid' or [] (default) - As described above.
% 
%           'same' - The region of interest covers the whole of V, and the
%           output vector has the same size as V. Reflection at the ends is
%           used to extrapolate V if it does not wrap around.
% 
%           A 2-element row vector with elements [MINI, MAXI] giving the
%           first and last indices of a section of V. The results vector
%           will have length MAXI-MINI+1 and will contain the smoothed
%           values for the specified section of the input vector.
%           Reflection at the ends is used for extrapolation of V when
%           necessary.
%
%       'Wrap' - Specifies whether the vector wraps round (i.e. has
%       periodic boundary conditions). The value must be a logical scalar,
%       or [] which is equivalent to false. If the vector wraps, the
%       default size of the output vector will be the same as the input
%       vector size. The default is false.
% 
%       'Centred' - Specifies whether centred differencing is used:
% 
%           true (default) - Gradient estimates are at element centres.
%           Thus VG(I) is a gradient estimate centred on V(I) if the 'same'
%           region specification is used. Otherwise it is centred on
%           V(I-REGION(1)+1). This is achieved by centred differencing,
%           using the two nearest neighbours.
% 
%           false - Gradient estimates are between adjacent elements of V.
%           Thus VG(I) is a gradient estimate centred on V(I-0.5) if the
%           'same' region specification is used. Otherwise it is centred on
%           V(I-REGION(1)+0.5). This is achieved by differencing adjacent
%           elements. The estimate is more localised than for centred
%           differencing, but is offset towards the origin relative to the
%           corresponding element of V.
%
% See also: gradients_n, gradients_xy, gradients_xyt, gsmooth

% Copyright David Young 2014

if ~isvector(x)
    error('DavidYoung:gradients_x:notVector', ...
        'First argument must be a vector');
end

row = isrow(x);
if row
    x = x.';
end

if length(varargin) > 1 && strcmp(varargin{2}, 'RegionOnly')
    % deal with region only case
    xg = gradients_n(x, varargin{:});
else
    [xg, region] = gradients_n(x, varargin{:});
    xg = xg{1};
    if row
        xg = xg.';
    end
end

end
