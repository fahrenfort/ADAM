function output = cond_string(varargin)
% create a comma separated string from the intersection of the input arrays
%
% By J.J.Fahrenfort, VU, 2015

output = varargin{1};
for c = 2:numel(varargin)
    output = intersect(output,varargin{c});
end
output = vec2str(output);