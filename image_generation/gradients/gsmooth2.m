function [a, reg] = gsmooth2(a, sigmas, varargin)
%GSMOOTH2 Gaussian image smoothing
%   [SMOOTH, REGION] = GSMOOTH2(A, SIGMAS) carries out Gaussian smoothing
%   on a 2-dimensional array. This differs from GSMOOTHN only in that image
%   coordinate ordering is used for anisotropic smoothing and boundary
%   conditions; the region of interest argument and result are, however, in
%   array coordinate order. A condensed argument form for backward
%   compatibility is allowed.
%  
%   A must be a matrix of class double or single.
%
%   SIGMAS specifies the smoothing constants. This may be a 2-vector of the
%   form [SIGMA_X, SIGMA_Y] giving the "sigma" parameters of the 1-D
%   Gaussian masks used for smoothing in the X and Y (row and column
%   respectively) directions. Alternatively SIGMAS may be a scalar, in
%   which case the same constant is used on both dimensions.
%
%       The output array will be smaller than the input array in order to
%       avoid having to extrapolate beyond the boundaries of the array. The
%       result REGION reports the region of the input array for which
%       corresponding output values have been computed, in the form
%       [MINROW, MAXROW, MINCOL, MAXCOL]. The size of the output array is
%       [MAXROW-MINROW+1, MAXCOL-MINCOL+1]. The reduction in size depends
%       on the smoothing parameters, and is chosen so that the smoothing
%       masks are a good approximation to the untruncated Gaussian mask.
% 
%   REGION = GSMOOTH2(A, SIGMAS, 'RegionOnly') returns only the region.
%
%   [SMOOTH, REGION] = GSMOOTH2(A, SIGMAS, NAME, VALUE, ...) allows
%   additional parameters to be set using name-value pairs. Names may be as
%   follows:
% 
%       'Region' - Allows a region of interest in A to be specified. Its
%       value may be:
% 
%           'valid' or [] - (default). As described above.
% 
%           'same' - The region of interest covers the whole of A, and the
%           output array has the same size as A. Reflection at the
%           boundaries is used to extrapolate A on dimensions where it does
%           not wrap around.
% 
%           A 4-element row vector with elements [MINROW, MAXROW, MINCOL,
%           MAXCOL] describing a rectangular volume. The results array will
%           have size [MAXROW-MINROW+1, MAXCOL-MINCOL+1] and will contain
%           the smoothed values for the specified region of the input
%           array. Reflection at the boundaries is used for extrapolation
%           of A when necessary.
% 
%       'Wrap' - Specifies that the array wraps round (i.e. has periodic
%       boundary conditions) on one or more axes. The value may be:
% 
%           A logical vector of length D of the form [WRAP_X, WRAP_Y]. If
%           WRAP_X is true the rows of A are circular; if WRAP_Y is true
%           the columns of A are circular. If a dimension wraps, the size
%           of that dimension will not be reduced in the default output
%           region.
% 
%           A logical scalar specifying that both dimensions are or are not
%           circular. The default is false.
% 
%           [] - Same as false.
% 
%   The following argument lists are allowed for backward compatibility:
% 
%   REGION = GSMOOTH2(A, SIGMAS, 'region') returns only the default region
% 
%   [SMOOTH, REGION] = GSMOOTH2(A, SIGMAS, REGION)
%       where REGION is empty, the string 'same', or a 4-element vector,
%       specifies a 'Region' parameter without the parameter name.
% 
%   [SMOOTH, REGION] = GSMOOTH2(A, SIGMAS, REGION, WRAP)
%       where WRAP is empty, a logical scalar, or a 2-element logical
%       vector, specifies a 'Wrap' parameter without the parameter name.
% 
% Example:
%
%       img = double(imread('pout.tif'));
%       imsmooth = gsmooth2(img, 4, 'Region', 'same');
% 
% See also: GSMOOTH, GSMOOTHN

% Copyright David Young 2014

[sigmas, reg, wraps, regonly] = checkinputs(a, sigmas, varargin{:});

% treat col as 2-D - gsmoothn would treat it as 1-D
col = iscolumn(a);
if col
    a = a.';
end

if regonly
    a = gsmoothn(a, sigmas, 'RegionOnly', 'Region', reg, 'Wrap', wraps);
    if col
        a = a([3 4 1 2]);
    end
else
    [a, reg] = gsmoothn(a, sigmas, 'Region', reg, 'Wrap', wraps);
    if col
        a = a.';
    end
end

end


function [sigmas, reg, wraps, regonly] = checkinputs(a, sigmas, varargin)
% Check arguments and get defaults, change convention. 
% Some checking will be done in gsmoothn, no need to repeat it.

% array argument
if ~ismatrix(a)
    error('DavidYoung:gsmooth2:notMatrix', ...
        'First argument must be a matrix');
end

% argument checking to allow switching between old and new conventions
% no need for region elements to be positive
checkreg = @(r) isempty(r) || ...
    checkattributes(r, {'numeric'}, {'integer' 'size' [1 4]});
checkwrap = @(w) isempty(w) || ...
    (islogical(w) && (isscalar(w) || isequal(size(w), [1 2])));

inp = inputParser;
% non param-val form for arguments, plus RegionOnly option
inp.addOptional('regnp', [], @(r) checkreg(r) || ...
    (ischar(r) && ismember(r, {'same', 'region', 'RegionOnly'})));
inp.addOptional('wrapnp', [], checkwrap);

% param-val form
inp.addParameter('Region', [], @(r) checkreg(r) || ...
    (ischar(r) && ismember(r, {'same', 'valid'})));
inp.addParameter('Wrap', [], checkwrap);

inp.parse(varargin{:});

% Check for consistency and select argument values
if ismember('wrapnp', inp.UsingDefaults) && ...
        (ismember('regnp', inp.UsingDefaults) || ...
        strcmp(inp.Results.regnp, 'RegionOnly'))
    % using name-val pairs
    regonly = ~isempty(inp.Results.regnp);
    reg = inp.Results.Region;
    wraps = inp.Results.Wrap;
elseif all(ismember({'Region' 'Wrap'}, inp.UsingDefaults))
    % using old arguments
    reg = inp.Results.regnp;
    regonly = strcmp(reg, 'region');
    if regonly
        reg = [];
    end
    wraps = inp.Results.wrapnp;
else
    error('DavidYoung:gsmooth2:inconsistentArgStyle', ...
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
