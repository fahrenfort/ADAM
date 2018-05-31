function runIntersect = LCS(varargin)
% LCS find longest common substring in a series of
% LCS(stringA,stringB,stringC,<...>);
% Using a cell array of strings you can use:
% LCS(cellstring{:});
% calls LCSSubStr (see below)
% 
if isempty(varargin),
    error('No inputs specified.')
else
    setArray = varargin;
end

runIntersect = setArray{1};
for i = 2:length(setArray),
    runIntersect = LCSubstr(runIntersect,setArray{i});
    if isempty(runIntersect),
        return
    end
end

function [ret,z,L]=LCSubstr(s,t)
%LCSSubStr is a function that return two string's largest common part
% Dr. XU Zhiping, School of Computer Scicence, Fudan University
% --------------------------------------------------------------------
%  
%	Usage:
%	 [ret]=LCSubstr(s,t)
%	s: input string 1
%	t: input string 2
%	ret: Largest Common String
%   z:   Largest Common String Length
%   L :  Compare matrix 
%   Example:
%   >>a='This is  very common string';
%   >>b='string is very common';
%   >>[ret]=LCSubstr(a,b)
%   ret =
%
%  very common

m=length(s);
n=length(t);
L=zeros(m+1,n+1);
z=1;
ret=[];
for i=1:m
    for j=1:n
        if s(i)==t(j)
            if i==1 || j==1
                L(i,j)=1;
            else
                L(i,j)=L(i-1,j-1)+1;
            end
            if L(i,j)>z
                z=L(i,j);
                ret=[];
            end
            if L(i,j)==z
                if isempty(ret)
                  ret=  s(i-z+1:i);                
                else
                  ret=[ret,',',s(i-z+1:i)];
                end
            end
        end
    end
end