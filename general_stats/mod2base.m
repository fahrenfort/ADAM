function outval = mod2base(inval,base)
% function to compute the numbers in inval, cut off to a maximum number
% that is defined in base. Inval can be a scalar or a vector, base must be
% a scalar. Examples:
% mod2base(34,30) returns 4
% mod2base([23,25],24) returns [23,1].
% mod2base([1002 10983], 0) returns [1002 10983]
% mod2base(48,24) returns 24
%
% By J.J.Fahrenfort, VU 2016

if nargin<2
    error('this function requires at least two arguments');
end
outval = mod(inval,base);
outval(outval == 0) = base;