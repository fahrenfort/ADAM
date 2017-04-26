function [xg, yg, region] = gradients_xy(im, varargin)
%GRADIENTS_XY Estimate gradients in 2 dimensions with Gaussian smoothing
%   [XG, YG, REGION] = GRADIENTS_XY(A) estimates the gradients of A using
%   centred differences.  This differs from GRADIENTS_N only in that image
%   coordinate ordering is used for anisotropic smoothing, boundary
%   conditions and the results arrays; the region of interest is, however,
%   in array coordinate order. A condensed argument form for backward
%   compatibility is allowed.
% 
%       A is a matrix of class double or single.
%
%       XG and YG are arrays of estimates of the gradients, computed using
%       symmetric differencing between the two nearest neighbours of each
%       voxel on the axis under consideration. XG is the gradient along the
%       rows of A, YG is the gradient along the columns of A. Each array
%       will be smaller than A to avoid needing to extrapolate A.
% 
%       REGION is the region of A for which gradients have been computed,
%       given as a 4-element vector with the minimum and maximum indices
%       for each dimension. In the non-smoothing case it will be [2
%       size(A,1)-1 2 size(A,2)-1].
% 
%   [XG, YG, REGION] = GRADIENTS_XY(A, SIGMAS) carries out Gaussian
%   smoothing and differencing to estimate the gradients of A along each
%   dimension.
%
%       SIGMAS specifies smoothing constants. It may be:
% 
%           A 2-vector of the form [SIGMAX SIGMAY]. SIGMAX and SIGMAY are
%           the "sigma" parameters of the 1D Gaussian masks used for
%           smoothing along the rows and columns of A respectively. (Note
%           that the order is different to that used in GRADIENTS_N.)
% 
%           A scalar, in which case the same constant is used on both
%           dimensions. If SIGMAS is zero, no smoothing is done.
%
%       The output arrays will be smaller than the input array in order to
%       avoid having to extrapolate beyond the boundaries of the array. The
%       result REGION reports the region of the input array for which the
%       output values have been estimated, in the form [MINROW, MAXROW,
%       MINCOL, MAXCOL]. The size of the output arrays is [MAXROW-MINROW+1,
%       MAXCOL-MINCOL+1]. The reduction in size depends on the smoothing
%       parameters, and is chosen so that the smoothing masks are a good
%       approximation to the untruncated Gaussian mask.
%
%   REGION = GRADIENTS_XY(A, SIGMAS, 'RegionOnly') returns only the region.
% 
%   [XG, YG, REGION] = GRADIENTS_XY(A, SIGMAS, NAME, VALUE, ...)
%   [XG, YG, REGION] = GRADIENTS_XY(A, NAME, VALUE, ...) 
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
%           A 4-element row vector with elements [MINROW, MAXROW, MINCOL,
%           MAXCOL] describing a rectangular region. The results arrays
%           will have size [MAXROW-MINROW+1, MAXCOL-MINCOL+1] and will
%           contain the smoothed values for the specified region of the
%           input array. Reflection at the boundaries is used for
%           extrapolation of A when necessary.
%
%       'Wrap' - Specifies that the array wraps round (i.e. has periodic
%       boundary conditions) on one or more axes. The value may be:
% 
%           A logical 2-vector of the form [WRAPX WRAPY]. If WRAPX is true
%           the rows of A are circular, if WRAPY is true the columns of A
%           are circular. If a dimension wraps, the size of that dimension
%           will not be reduced in the default output region. (Note that
%           the ordering is different to GRADIENTS_N.)
% 
%           A logical scalar specifying that all dimensions are or are not
%           circular. The default is false.
% 
%           [] - Same as false.
% 
%       'Centred' - Specifies whether centred differencing is used:
% 
%           true (default) - Gradient estimates are at pixel centres. Thus
%           XG(I,J) and YG(I,J) are gradient estimates centred on A(I,J) if
%           the 'same' region specification is used. Otherwise they are
%           centred on A(I-REGION(1)+1, J-REGION(3)+1). This is achieved by
%           centred differencing, using the two nearest neighbours along
%           the row or column dimension.
% 
%           false - Gradient estimates are at pixel corners. Thus XG(I,J)
%           and YG(I,J) are gradient estimates centred on A(I-0.5, J-0.5)
%           if the 'same' region specification is used. Otherwise they are
%           centred on A(I-REGION(1)+0.5, J-REGION(3)+0.5). This is
%           achieved by differencing the average values on opposite faces
%           of a pixel-sized square. The estimate is more localised than
%           for centred differencing, but is offset towards the origin
%           relative to the corresponding element of A.
%
%   The following argument lists are allowed for backward compatibility:
% 
%   REGION = GRADIENTS_XY(A, SIGMAS, 'region') returns only the default
%   region
% 
%   [XG, YG, REGION] = GRADIENTS_XY(A, SIGMAS, REGION)
%       where REGION is empty, the string 'same', or a 4-element vector,
%       specifies a 'Region' parameter without the parameter name.
% 
%   [XG, YG, REGION] = GRADIENTS_XY(A, SIGMAS, REGION, WRAP)
%       where WRAP is empty, a logical scalar, or a 2-element logical
%       vector, specifies a 'Wrap' parameter without the parameter name.
% 
% Example:
%     
%         img = double(imread('pout.tif'));
%         [xg, yg] = gradients_xy(img, 4, 'Region', 'same');
% 
% See also: gradients_x, gradients_n, gradients_xyt, gsmooth2

% Copyright David Young 2014

% check arguments, change convention
[sigmas, region, wraps, regonly, centred] = checkinputs(im, varargin{:});

% treat col as 2-D - gradients_n would treat it as 1-D
col = iscolumn(im);
if col
    im = im.';
end

if regonly
    xg = gradients_n(im, sigmas, 'RegionOnly', ...
        'Region', region, 'Wrap', wraps, 'Centred', centred);
    if col
        xg = xg([3 4 1 2]);
    end
else
    [g, region] = gradients_n(im, sigmas, ...
        'Region', region, 'Wrap', wraps, 'Centred', centred);
    xg = g{2};
    yg = g{1};
    if col
        xg = xg.';
        yg = yg.';
    end
end

end

% -------------------------------------------------------------------------

function [sigmas, reg, wraps, regonly, centred] = checkinputs(a, varargin)
% Check arguments and get defaults, change convention. 
% Some checking will be done in gradients_n, no need to repeat it.

% array argument
if ~ismatrix(a)
    error('DavidYoung:gradients_xy:notMatrix', ...
        'First argument must be a matrix');
end

% argument checking to allow switching between old and new conventions
% no need for region elements to be positive
checkreg = @(r) isempty(r) || ...
    checkattributes(r, {'numeric'}, {'integer' 'size' [1 4]});
checkwrap = @(w) isempty(w) || ...
    (islogical(w) && (isscalar(w) || isequal(size(w), [1 2])));

inp = inputParser;
inp.addOptional('sigma', 0);
% non param-val form for arguments, plus RegionOnly option
inp.addOptional('regnp', [], @(r) checkreg(r) || ...
    (ischar(r) && ismember(r, {'same', 'region', 'RegionOnly'})));
inp.addOptional('wrapnp', [], checkwrap);

% param-val form
inp.addParameter('Region', [], @(r) checkreg(r) || ...
    (ischar(r) && ismember(r, {'same', 'valid'})));
inp.addParameter('Wrap', [], checkwrap);
inp.addParameter('Centred', true);

inp.parse(varargin{:});

sigmas = inp.Results.sigma;
centred = inp.Results.Centred;

% Check for consistency and select argument values
if ismember('wrapnp', inp.UsingDefaults) && ...
        (ismember('regnp', inp.UsingDefaults) || ...
        strcmp(inp.Results.regnp, 'RegionOnly'))
    % using name-val pairs
    regonly = ~isempty(inp.Results.regnp);
    reg = inp.Results.Region;
    wraps = inp.Results.Wrap;
elseif all(ismember({'Region' 'Wrap' 'Centred'}, inp.UsingDefaults))
    % using old arguments
    reg = inp.Results.regnp;
    regonly = strcmp(reg, 'region');
    if regonly
        reg = [];
    end
    wraps = inp.Results.wrapnp;
else
    error('DavidYoung:gradients_xy:inconsistentArgStyle', ...
        'Inconsistent use of name-value pairs and other arguments');
end

% switch sigmas and wraps from image to array convention
if isequal(numel(sigmas), 2)
    sigmas = sigmas([2 1]);
end
if isequal(numel(wraps), 2)
    wraps = wraps([2 1]);
end

end