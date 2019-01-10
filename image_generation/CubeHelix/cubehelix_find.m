function [params,irange,domain,resnorm] = cubehelix_find(map,posn,prms)
% Return parameter values for the CubeHelix colorscheme that matches the input colormap.
%
% (c) 2013 Stephen Cobeldick
%
% So you have a nice CubeHelix colormap but you can't remember the exact
% parameter values that were used to define it, or perhaps you want to match
% a CubeHelix scheme to your document's chosen colorscheme: this function can
% help! An optimization routine finds the CubeHelix colorscheme that best matches
% the input RGB colormap, and returns its parameter values (takes a few seconds).
%
%%% Syntax:
%  [params,irange,domain] = cubehelix_find(map)
%  [params,irange,domain] = cubehelix_find(map,posn)
%  [params,irange,domain] = cubehelix_find(map,posn,prms)
%
% Note1: Requires the Optimization Toolbox function LSQNONLIN.
% Note2: The parameter <start> is modulus three, i.e. 3==0.
% Note3: The search bounds are lb=[0,-20,0,0,0,0,0,0] and ub=[3,+20,3,3,1,1,1,1];
%
% See also CUBEHELIX BREWERMAP RGBPLOT COLORMAP LSQNONLIN OPTIMSET
%
%% Matching a Colorscheme %%
%
% The basic syntax assumes that the input colormap is a CubeHelix colormap:
%  [params,irange,domain] = cubehelix_find(map) % map = cubehelix(...)
%
% You can also find the closest CubeHelix colorscheme that matches a
% colormap of arbitrary colors (e.g. a corporate or document colorscheme).
% This can be achieved by supplying a colormap of the desired colors, and a
% vector giving their relative positions in the final CubeHelix colorscheme:
%
% map  = [0.16,0.59,0.87; 0.78,0.55,0.35; 0.78,0.9,0.9]; % three nodes
% posn = [0.25,0.5,0.75]; % relative positions of the three <map> nodes
% [params,irange,domain] = cubehelix_find(map,posn);
% cubehelix_view(5,params,irange,domain)
%
%% Holding Parameter Values Constant %%
%
% In some situations it may be desirable to force some parameter values to
% stay constant. The third (optional) input can be used to specify any of
% the eight parameter values, other values to be determined are set to NaN:
%
% prms = [NaN,NaN,NaN,1,NaN,NaN,0,1]; % constant gamma=1 and domain=[0,1]
% [params,irange,domain] = cubehelix_find(map,posn,prms);
% cubehelix_view(5,params,irange,domain)
%
%% Examples %%
%
% cubehelix_find(cubehelix(10))
%  ans = [0.5,-1.5,1,1]
%
% map = cubehelix(10, 1.4,-0.7,0.9,1.2, [0.05,0.97]);
% [params,irange,domain] = cubehelix_find(map)
%  params = [1.4,-0.7,0.9,1.2]
%  irange = [0.05,0.97]
%  domain = [0,1]
%
% map = cubehelix(64, [2.3,0.4,0.5,0.6], [0.05,0.24], [0.19,0.85]);
% [params,irange,domain] = cubehelix_find(map)
%  params = [2.3,0.4,0.5,0.6]
%  irange = [0.05,0.24]
%  domain = [0.19,0.85]
%
%% Input and Output Arguments %%
%
%%% Inputs (*=default):
%  map  = NumericMatrix, an RGB colormap to be matched to a CubeHelix colorscheme.
%  posn = NumericVector, the relative positions of the nodes in <map>. Size 1xrows(map).
%  prms = NumericVector, vector of CubeHelix parameters to keep constant:
%         =[start,rots,sat,gamma,irange,domain]. NaN indicates non-constant values.
%%% Outputs:
%  params = NumericVector, the CubeHelix parameters derived from <map>: [start,rots,sat,gamma].
%  irange = NumericVector, range of brightness levels of the scheme's endnodes. Size 1x2.
%  domain = NumericVector, domain of the cubehelix calculation (endnode positions). Size 1x2.
%
% [params,irange,domain] = cubehelix_find(map,*posn,*prms)

%% Input Wrangling %%
%
ism = @(x)all(0<=x(:))&&all(x(:)<=1);
isn = @(x,n)isnumeric(x)&&isreal(x);
%
assert(isn(map)&&ismatrix(map)&&size(map,2)==3&&ism(map),...
	'Input argument <map> must be a colormap of RGB values (size Nx3).')
% Grayscale colormap equivalent:
map = double(map);
mag = map*[0.298936;0.587043;0.114021];
%
R = size(map,1);
%
if nargin>2
	assert(isn(prms)&&isvector(prms)&&all(~isinf(prms(:)))&&numel(prms)==8,...
		'Input argument <prms> must be a real numeric vector of size 1x8.')
	idx = isnan(prms);
	if ~any(idx)
		params = prms(1:4);
		irange = prms(5:6);
		domain = prms(7:8);
		resnorm = NaN;
		return
	end
	X = @(p)p(idx);
else
	idx = true(1,8);
	X = @(p)p;
end
%
%% LSQNONLIN Parameter Bounds %%
%
if R>2 % default algorithm cannot solve <3 colormap rows
	lb = X([0,-20,0,0,0,0,0,0]);
	ub = X([3,+20,3,3,1,1,1,1]);
	opt = optimset('Display','off', 'Algorithm','trust-region-reflective');
else % no bounds allowed for this algorithm
	lb = [];
	ub = [];
	opt = optimset('Display','off', 'Algorithm','levenberg-marquardt');
end
%
%% Indexing for Sample Colors %%
%
if nargin>1 && ~isempty(posn) % sampled colors
	assert(isn(posn)&&isvector(posn)&&ism(posn)&&numel(posn)==R,...
		'Input argument <posn> is a vector of the relative position in the colorscheme of each <map> node.')
	[N,D] = rat(posn);
	P = prod(D);
	Y = N.*(P./D);
	E = @(m)m(Y,:);
else % complete CubeHelix colormap
	P = R;
	E = @(m)m;
end
%
%% Define Solver Function %%
%
G = @(p)map-E(cubehelix(P,p(1:4),p(5:6),p(7:8)));
if nargin>2
	F = @(p)G(chfSubs(real(p),prms,idx));
else
	F = @(p)G(real(p));
end
%
%% Find Optimum %%
%
% Estimate rotations from the number of peaks&troughs:
D = diff(bsxfun(@minus,map,mag),1,1);
idy = D(1:end-1,:)>=0 & D(2:end,:)<0;
idz = D(1:end-1,:)<0 & D(2:end,:)>=0;
rots = ceil(mean(sum([idy,idz],1)));
%
% Range of start and rotation values to try:
Sm = linspace(0,3-eps,(1+rots)*7);
Rm = rots*[1;-1]*[0.7,1,1.3];
Rm = [0;Rm(:)];
[Sm,Rm] = meshgrid(Sm(1:max(1,idx(1)*end)),Rm(1:max(1,idx(2)*end)));
% Solve!
[Zc,rn] = arrayfun(@(s,r)lsqnonlin(F,X([s,r,1,1,0,1,0,1]),lb,ub,opt), Sm,Rm, 'UniformOutput',false);
% Pick the best solution:
[resnorm,idz] = min([rn{:}]);
if nargin>2
	out(idx) = real(Zc{idz});
	out(~idx) = prms(~idx);
else
	out = real(Zc{idz});
end
domain = out(7:8);
irange = out(5:6);
params = out(1:4);
params(1) = mod(params(1),3);
% Residual norm per row:
resnorm = resnorm / R;
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cubehelix_find
function q = chfSubs(p,v,idx)
% Zip together constant and non-constant parameter value vectors.
q(idx) = p;
q(~idx) = v(~idx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%chfSubs
% Copyright (c) 2013 Stephen Cobeldick
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license