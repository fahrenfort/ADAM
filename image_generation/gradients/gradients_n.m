function [g, region] = gradients_n(a, varargin)
%GRADIENTS_N Estimate gradients in N dimensions with Gaussian smoothing
%   [G, REGION] = GRADIENTS_N(A) estimates the gradients of A using centred
%   differences.
% 
%       A is an N-D array of class double or single. Its dimensionality D
%       is given by NDIMS(A) unless A is a column vector, when D = 1.
%
%       G is a cell array of estimates of the gradients, computed using
%       symmetric differencing between the two nearest neighbours of each
%       voxel on the axis under consideration. G{1} is the gradient along
%       the first dimension of A, etc. Each array will be smaller than A to
%       avoid needing to extrapolate A.
% 
%       REGION is the region of A for which gradients have been computed,
%       given as a 2*D vector with the minimum and maximum indices for each
%       dimension. In the non-smoothing case it will be [2 size(A,1)-1 2
%       size(A,2)-1 ...].
% 
%   [G, REGION] = GRADIENTS_N(A, SIGMAS) carries out Gaussian smoothing and
%   differencing to estimate the gradients of A along each dimension.
%
%       SIGMAS specifies smoothing constants. It may be:
% 
%           A vector of length D of the form [SIGMA1, SIGMA2, ...]. SIGMA1,
%           SIGMA2 etc. are the "sigma" parameters of the 1D Gaussian masks
%           used for smoothing along each dimension of A. (Note that the
%           order of the first two constants is different to that used in
%           GRADIENTS_XY and GRADIENTS_XYT.)
% 
%           A scalar, in which case the same constant is used on all
%           dimensions. If SIGMAS is zero, no smoothing is done.
%
%       The output arrays will be smaller than the input array in order to
%       avoid having to extrapolate beyond the boundaries of the array. The
%       result REGION reports the region of the input array for which the
%       output values have been estimated, in the form [MIN1, MAX1, MIN2,
%       MAX2, ...] where MIN1 and MAX1 are the limits of the region on the
%       first dimension, etc. The size of the output arrays is
%       [MAX1-MIN1+1, MAX2-MIN2+1, ...]. The reduction in size depends on
%       the smoothing parameters, and is chosen so that the smoothing masks
%       are a good approximation to the untruncated Gaussian mask.
%
%   REGION = GRADIENTS_N(A, SIGMAS, 'RegionOnly') returns only the region.
% 
%   [G, REGION] = GRADIENTS_N(A, SIGMAS, NAME, VALUE, ...)
%   [G, REGION] = GRADIENTS_N(A, NAME, VALUE, ...) 
%   set additional parameters using name-value pairs. Names may be:
% 
%       'Region' - Specifies a region of interest in A. Its value may be:
% 
%           'valid' or [] (default) - As described above.
% 
%           'same' - The region of interest covers the whole of A, and the
%           output array has the same size as A. Reflection at the
%           boundaries is used to extrapolate A on dimensions where it does
%           not wrap around.
% 
%           A 2*D-element row vector with elements [MIN1, MAX1, MIN2, MAX2,
%           ...] describing a rectangular volume. The results arrays will
%           have size [MAX1-MIN1+1, MAX2-MIN2+1, ...] and will contain the
%           smoothed values for the specified region of the input array.
%           Reflection at the boundaries is used for extrapolation of A
%           when necessary.
%
%       'Wrap' - Specifies that the array wraps round (i.e. has periodic
%       boundary conditions) on one or more axes. The value may be:
% 
%           A logical vector of length D of the form [WRAP1, WRAP2, ...].
%           If WRAP1 is true the first dimension is circular, etc. If a
%           dimension wraps, the size of that dimension will not be reduced
%           in the default output region.
% 
%           A logical scalar specifying that all dimensions are or are not
%           circular. The default is false.
% 
%           [] - Same as false.
% 
%       'Centred' - Specifies whether centred differencing is used:
% 
%           true (default) - Gradient estimates are at voxel centres. Thus
%           G{Q}(I,J,...) is a gradient estimate centred on A(I,J,...) if
%           the 'same' region specification is used. Otherwise it is
%           centred on A(I-REGION(1)+1, J-REGION(3)+1, ...). This is
%           achieved by centred differencing, using the two nearest
%           neighbours along the Q dimension.
% 
%           false - Gradient estimates are at voxel corners. Thus
%           G{Q}(I,J,...) is a gradient estimate centred on A(I-0.5, J-0.5,
%           ...) if the 'same' region specification is used. Otherwise it
%           is centred on A(I-REGION(1)+0.5, J-REGION(3)+0.5, ...). This is
%           achieved by differencing the average values on opposite faces
%           of a voxel-sized hypercube. The estimate is more localised than
%           for centred differencing, but is offset towards the origin
%           relative to the corresponding element of A.
%
% See also: gradients_x, gradients_xy, gradients_xyt, gsmoothn

% Copyright David Young 2014

% check arguments and get defaults
[sigmas, wraps, region, regonly, d, symdiff, dmargin] = ...
    checkinputs(a, varargin{:});

% can stop now if only the region to be returned
if regonly
    g = region;
    return;
end
    
% expand the region to allow for subsequent differencing operation
regdiff = region + repmat(dmargin, 1, d);

% region selection and spatial smoothing (do this before sum and difference
% to reduce processing if region is small)
asmth = gsmoothn(a, sigmas, 'Region', regdiff, 'Wrap', wraps);

% Differencing mask
if symdiff
    % get the average gradient over 2 voxels so centred
    dmask = [1; 0; -1]/2;
    % Index ranges for trimming undifferenced dimensions
    sz = size(asmth);
    sz = sz(1:d);     % omit trailing 1 in 1-D case
    triminds = arrayfun(@(x) 2:x-1, sz, 'UniformOutput', false);
else
    % take adjacent voxels, not centred on voxel
    dmask = repmat([1; -1]/(2^(d-1)), [1 2+zeros(1,d-1)]);
end

g = cell(1, d);
perm = (1:max(2,d)).';   % permute needs at least 2 dimensions to order

for i = 1:d
    
    % orient mask along current dimension
    dmaskr = permute(dmask, circshift(perm, i-1));
    
    if symdiff
        % trim other dimensions so outputs all same size
        t = triminds;
        t{i} = 1:sz(i);        %
        g{i} = convn(asmth(t{:}), dmaskr, 'valid');
    else
        g{i} = convn(asmth, dmaskr, 'valid');
    end
    
    
end

end

% -------------------------------------------------------------------------

function [sigmas, wraps, region, regonly, d, symdiff, dmargin] = ...
    checkinputs(a, varargin)
% Check arguments and get defaults

% Most checking done in gsmoothn, no need here
inp = inputParser;
inp.addOptional('sigma', 0);
inp.addOptional('regonly', '', @(s) strcmp(s, 'RegionOnly'));
inp.addParameter('Region', []);
inp.addParameter('Wrap', []);
inp.addParameter('Centred', true);
inp.parse(varargin{:});
regonly = ~isempty(inp.Results.regonly);
sigmas = inp.Results.sigma;
region = inp.Results.Region;
wraps = inp.Results.Wrap;
symdiff = inp.Results.Centred;

% number of dimensions, defined as number of last non-singleton dimension
if iscolumn(a)
    d = 1;
else
    d = ndims(a);
end

if isempty(wraps)
    wraps = false(1, d);
elseif isscalar(wraps)
    wraps = repmat(wraps, 1, d);
end

if symdiff
    dmargin = [-1 1];   % region margin for differencing
else
    dmargin = [-1 0];
end

if isempty(region) || strcmp(region, 'valid')
    % default region - small enough not to need extrapolation
    region = gsmoothn(a, sigmas, 'RegionOnly', 'Wrap', wraps);
    % -dmargin because contracting to get output region
    region = region + ...
        repmat(-dmargin, 1, d) .* double(~reshape([wraps; wraps], 1, 2*d));
elseif strcmp(region, 'same')
    sz = size(a);
    region = reshape([ones(1, d); sz(1:d)], 1, 2*d);
end
if any(region(2:2:end) < region(1:2:end))
    error('DavidYoung:gradients_n:badreg', ...
        'REGION or array size too small');
end

end

