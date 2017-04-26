function [value,row,colomn]=max2d(x)
% [value,row,colomn]=min2d(x)
% reutrns tha value of the first maximum element in a 2-D matrix 'x' and 
% its position (row & colomn)
%
% By: Abdulrahman Ikram Siddiq
% Kirkuk - IRAQ
% Wednsday Nov.9th 2011 10:23 PM

[w,j]=max(x);
[value,i]=max(w);
colomn=i;
row=j(i);