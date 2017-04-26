function [value,row,colomn]=min2d(x)
% [value,row,colomn]=min2d(x)
% reutrns tha value of the first minimum element in a 2-D matrix 'x' and 
% its position (row & colomn)
%
% By: Abdulrahman Ikram Siddiq
% Kirkuk - IRAQ
% Wednsday Nov.9th 2011 10:20 PM
[w,j]=min(x);
[value,i]=min(w);
colomn=i;
row=j(i);