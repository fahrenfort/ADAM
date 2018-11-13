function [C,IA,IC,N] = nt_unique(A, varargin)
%[C,IA,IC,N] = nt_unique(A, varargin) - unique with counts
%
%  N: number of occurrences
%
% See unique for definition of other variables.
% 
% NoiseTools

[Asorted,iSort]=sortrows(A);
[~,iReverse]=sort(iSort);
[C,IA,IC]=unique(Asorted,varargin{:});
p=find([1;diff(IC);1]);  
N=diff(p); 

IA=iSort(IA);
IC=IC(iReverse);