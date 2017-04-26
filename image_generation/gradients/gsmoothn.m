function [a, reg] = gsmoothn(a, sigmas, varargin)
%GSMOOTHN Gaussian N-dimensional array smoothing
%   [SMOOTH, REGION] = GSMOOTHN(A, SIGMAS) carries out Gaussian
%   smoothing on an N-dimensional array.
% 
%   A must be an array of class double or single. Its dimensionality D is
%   given by NDIMS(A) unless A is a column vector, when D = 1.
%
%   SIGMAS specifies smoothing constants. This may be:
% 
%       A vector of length D of the form [SIGMA1, SIGMA2, ...]. SIGMA1,
%       SIGMA2 etc. are the "sigma" parameters of the 1D Gaussian masks
%       used for smoothing along each dimension of A. (Note that the order
%       of the first two constants is different to that used in
%       GRADIENTS_XY and GRADIENTS_XYT.)
% 
%       A scalar, in which case the same constant is used on all
%       dimensions.
%
%   SMOOTH is the output array of smoothed values.
% 
%       The output array will be smaller than the input array in order to
%       avoid having to extrapolate beyond the boundaries of the array. The
%       result REGION reports the region of the input array for which
%       corresponding output values have been computed, in the form [MIN1,
%       MAX1, MIN2, MAX2, ...] where MIN1 and MAX1 are the limits of the
%       region on the first dimension, etc. The size of the output array is
%       [MAX1-MIN1+1, MAX2-MIN2+1, ...]. The reduction in size depends on
%       the smoothing parameters, and is chosen so that the smoothing masks
%       are a good approximation to the untruncated Gaussian mask.
% 
%   REGION = GSMOOTHN(A, SIGMAS, 'RegionOnly') returns only the region.
%
%   [SMOOTH, REGION] = GSMOOTHN(A, SIGMAS, NAME, VALUE, ...) sets
%   additional parameters using name-value pairs. Names may be:
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
%           ...] describing a rectangular volume. The results array will
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
% See also: GSMOOTH, GSMOOTH2

% Copyright David Young 2014

[sigmas, bcons, reg, convreg, ronly, d] = ...
    checkinputs(a, sigmas, varargin{:});

if ronly
    a = reg;
else
    if ~isempty(convreg)            % trim/pad array
        lims = mat2cell(convreg, 1, repmat(2, 1, d));
        lims = cellfun(@(x) x(1):x(2), lims, 'UniformOutput', false);
        indargs = cell(1, 2*d);
        indargs(1:2:end) = lims;
        indargs(2:2:end) = bcons;
        
        a = exindex(a, indargs{:});
    end
    for i = 1:d
        a = gsmooth1(a, i, sigmas(i));
    end
end

end

% -------------------------------------------------------------------------

function [sigmas, bcons, reg, convreg, regonly, d] = ...
    checkinputs(a, sigmas, varargin)
% Check arguments and get defaults, plus input/output convolution regions
% and boundary conditions

% array argument
validateattributes(a, {'double' 'single'}, {}, mfilename, 'a');

% number of dimensions, defined as number of last non-singleton dimension
if iscolumn(a)
    d = 1;
else
    d = ndims(a);
end

% sigmas argument
validateattributes(sigmas, {'double'}, {'nonnegative' 'real' 'vector'}, ...
    mfilename, 'sigmas');
if isscalar(sigmas)
    sigmas = repmat(sigmas, 1, d);
elseif ~isequal(length(sigmas), d)
    error('DavidYoung:gsmoothn:badsigmas', 'SIGMAS wrong size');
end

inp = inputParser;
inp.addOptional('regonly', '', @(s) strcmp(s, 'RegionOnly'));
% no need for region elements to be positive
inp.addParameter('Region', [], @(r) ...
    isempty(r) || ...
    (ischar(r) && ismember(r, {'valid' 'same'})) || ...
    checkattributes(r, {'numeric'}, {'integer' 'size' [1 2*d]}));
inp.addParameter('Wrap', [], @(w) ...
    isempty(w) || ...
    (islogical(w) && (isscalar(w) || isequal(size(w), [1 d]))));
inp.parse(varargin{:});

% wraps argument
wraps = inp.Results.Wrap;
if isempty(wraps)
    wraps = false(1, d);
elseif isscalar(wraps)
    wraps = repmat(wraps, 1, d);
end
boundopts = {'symmetric' 'circular'};
bcons = boundopts(wraps+1);

% region argument
regonly = ~isempty(inp.Results.regonly);
reg = inp.Results.Region;

% whole array region
sz = size(a);
imreg = reshape([ones(1, d); sz(1:d)], 1, 2*d);

% convolution margins and wrap multipliers
mrg = gausshsize(sigmas);
mrg = reshape([mrg; -1*mrg], 1, 2*d);

if isempty(reg) || strcmp(reg, 'valid')
    % default region - small enough not to need extrapolation - shrink on
    % non-wrapped dimensions
    reg = imreg + mrg .* double(~reshape([wraps; wraps], 1, 2*d));
elseif strcmp(reg, 'same')
    reg = imreg;
end
if any(reg(2:2:end) < reg(1:2:end))
    error('DavidYoung:gsmoothn:badreg', 'REGION or array size too small');
end
% compute input region for convolution - expand on all dimensions
convreg = reg - mrg;    % expand
if isequal(convreg, imreg)
    convreg = [];   % signal no trimming or padding
end

end

% -------------------------------------------------------------------------

function a = gsmooth1(a, dim, sigma)
% Smooth an array A along dimension DIM with a 1D Gaussian mask of
% parameter SIGMA

if sigma > 0
    mlen = 2*gausshsize(sigma) + 1;  % reasonable truncation    
    mask = fspecial('gauss', [mlen 1], sigma);
    msize = ones(1, ndims(a));
    msize(dim) = mlen;
    mask = reshape(mask, msize);
    
    a = convn(a, mask, 'valid');
end

end

% -------------------------------------------------------------------------

function hsize = gausshsize(sigma)
% Default for the limit on a Gaussian mask of parameter sigma.
% Produces a reasonable degree of truncation without too much error.
hsize = ceil(2.6*sigma);
end

