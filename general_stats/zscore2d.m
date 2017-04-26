function normalized = zscore2d(averages)
% function to normalize data across a matrix
% e.g. takes an array that has the form data(subjectnr,conditionnr)
% that is, the mean and sd is taken across the entire matrix, but the
% columns and rows are retained in the output
%
% J.J.Fahrenfort, VU 2016
ntot = numel(averages);
origsize = size(averages);
normalized = reshape(zscore(reshape(averages,ntot,1)),origsize);