function [xg, yg, tg, region] = gradients_xyt(image1, image2, varargin)
%GRADIENTS_XYT Estimate spatial and temporal grey-level gradients
%   [XG, YG, TG, REGION] = GRADIENTS_XYT(IM1, IM2) estimates the spatial
%   gradients of the average of IM1 and IM2 using centred differences, and
%   the temporal gradient using the simple difference between the images.
% 
%       IM1 and IM2 are matrices that have the same size and are of class
%       double or single.
%
%       XG and YG are arrays of estimates of the spatial gradients,
%       computed using symmetric differencing between the two nearest
%       neighbours of each voxel on the axis under consideration, applied
%       to the average of IM1 and IM2. XG is the gradient along the rows,
%       YG is the gradient along the columns. Each array will be smaller
%       than IM1 to avoid needing to extrapolate.
% 
%       TG is an array of estimates of the temporal gradients, computed
%       simply as IM2-IM1, trimmed to the same size as XG and YG.
% 
%       REGION is the region of the input arrays for which gradients have
%       been computed, given as a 4-element vector with the minimum and
%       maximum indices for each dimension. In the non-smoothing case it
%       will be [2 size(IM1,1)-1 2 size(IM2,2)-1].
% 
%   [XG, YG, TG, REGION] = GRADIENTS_XYT(IM1, IM2, SIGMAS) carries out
%   Gaussian smoothing on the spatial dimensions prior to differencing.
%
%       SIGMAS specifies smoothing constants. It may be:
% 
%           A 2-vector of the form [SIGMAX SIGMAY]. SIGMAX and SIGMAY are
%           the "sigma" parameters of the 1D Gaussian masks used for
%           smoothing along the rows and columns respectively. (Note
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
%   REGION = GRADIENTS_XYT(IM1, IM2, SIGMAS, 'RegionOnly') returns only the
%   region.
% 
%   [XG, YG, TG, REGION] = GRADIENTS_XYT(IM1, IM2, SIGMAS, NAME, VALUE, ...)
%   [XG, YG, TG, REGION] = GRADIENTS_XYT(IM1, IM2, NAME, VALUE, ...) 
%   set additional parameters using name-value pairs. Names may be:
% 
%       'Region' - Specifies a spatial region of interest in the input
%       arrays. Its value may be:
% 
%           'valid' or [] (default) - As described above.
% 
%           'same' - The region of interest covers the whole of IM1, and
%           the output array has the same size as IM1. Reflection at the
%           boundaries is used to extrapolate the input arrays on
%           dimensions which do not wrap around.
% 
%           A 4-element row vector with elements [MINROW, MAXROW, MINCOL,
%           MAXCOL] describing a rectangular region. The results arrays
%           will have size [MAXROW-MINROW+1, MAXCOL-MINCOL+1] and will
%           contain the smoothed values for the specified region of the
%           input array. Reflection at the boundaries is used for spatial
%           extrapolation of IM1 and IM2 when necessary.
%
%       'Wrap' - Specifies that the arrays wrap round (i.e. have periodic
%       boundary conditions) on one or more axes. The value may be:
% 
%           A logical 2-vector of the form [WRAPX WRAPY]. If WRAPX is true
%           the rows are circular, if WRAPY is true the columns are
%           circular. If a dimension wraps, the size of that dimension will
%           not be reduced in the default output region. (Note that the
%           ordering is different to GRADIENTS_N.)
% 
%           A logical scalar specifying that all dimensions are or are not
%           circular. The default is false.
% 
%           [] - Same as false.
% 
%       'Centred' - Specifies whether centred differencing is used for the
%       spatial gradients:
% 
%           true (default) - Gradient estimates are at pixel centres. Thus
%           XG(I,J) and YG(I,J) are gradient estimates centred on IM1(I,J)
%           if the 'same' region specification is used. Otherwise they are
%           centred on IM1(I-REGION(1)+1, J-REGION(3)+1). This is achieved
%           by centred differencing, using the two nearest neighbours along
%           the row or column dimension.
% 
%           false - Gradient estimates are at pixel corners. Thus XG(I,J)
%           and YG(I,J) are gradient estimates centred on IM1(I-0.5, J-0.5)
%           if the 'same' region specification is used. Otherwise they are
%           centred on IM1(I-REGION(1)+0.5, J-REGION(3)+0.5). This is
%           achieved by differencing the average values on opposite faces
%           of a pixel-sized square. The estimate is more localised than
%           for centred differencing, but is offset towards the origin
%           relative to the corresponding element of IM1 or IM2.
% 
%           The temporal gradients are estimated by simple differencing,
%           and are always at the pixel centres.
%
%   The following argument lists are allowed for backward compatibility:
% 
%   REGION = GRADIENTS_XYT(IM1, IM2, SIGMAS, 'region') returns only the
%   default region
% 
%   [XG, YG, TG, REGION] = GRADIENTS_XYT(IM1, IM2, SIGMAS, REGION)
%       where REGION is empty, the string 'same', or a 4-element vector,
%       specifies a 'Region' parameter without the parameter name.
% 
%   [XG, YG, TG, REGION] = GRADIENTS_XYT(IM1, IM2, SIGMAS, REGION, WRAP)
%       where WRAP is empty, a logical scalar, or a 2-element logical
%       vector, specifies a 'Wrap' parameter without the parameter name.
% 
% Example:
%     
%   im1 = double(rgb2gray(imread('office_1.jpg')))/256;
%   im2 = double(rgb2gray(imread('office_2.jpg')))/256;
%   [xg, yg, tg] = gradients_xyt(im1, im2, 4, 'Region', 'same');
% 
% See also: gradients_x, gradients_n, gradients_xy, gsmooth2

% Copyright David Young 2014

% check arguments and get defaults
[sigmas, region, wraps, regonly, centred] = ...
    checkinputs(image1, image2, varargin{:});

% expand region
region = gradients_xy(image1, sigmas, 'RegionOnly', ...
    'Region', region, 'Wrap', wraps, 'Centred', centred);

% can stop now if only the region to be returned
if regonly
    xg = region;
    return;
end
    
% expand the region to allow for subsequent differencing operation
if centred
    regdiff = region + [-1 1 -1 1];
else
    regdiff = region + [-1 0 -1 0];
end

% region selection and spatial smoothing (do this before sum and difference
% to reduce processing if region is small)
imsmth1 = gsmooth2(image1, sigmas, 'Region', regdiff, 'Wrap', wraps);
imsmth2 = gsmooth2(image2, sigmas, 'Region', regdiff, 'Wrap', wraps);

% Spatial gradient on mean image
imav = (imsmth1 + imsmth2)/2;
[xg, yg] = gradients_xy(imav, 'Centred', centred);

% temporal differencing
tg = imsmth2 - imsmth1;
if centred
    tg = tg(2:end-1, 2:end-1);
else
    tg = tg(2:end, 2:end);
end

end

% -------------------------------------------------------------------------

function [sigmas, reg, wraps, regonly, centred] = ...
    checkinputs(im1, im2, varargin)
% Check arguments and get defaults, change convention. 
% Some checking will be done in gradients_n, no need to repeat it.

% array arguments
if ~isequal(size(im1), size(im2))
    error('DavidYoung:gradients_xyt:badsize', ...
        'Size of image2 does not match image1');
end
if ~ismatrix(im1)
    error('DavidYoung:gradients_xyt:notMatrix', ...
        'First two arguments must be matrices');
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
    error('DavidYoung:gradients_xyt:inconsistentArgStyle', ...
        'Inconsistent use of name-value pairs and other arguments');
end

end

