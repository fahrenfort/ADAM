function [y]=nt_mat2trial(x,w) 
%[y]=nt_trial2mat(x) - convert 3D matrix to trial cell array
%
%  y: trial array (each trial is channels * samples)
%
%  x: matrix (samples * channels * trials)
%  w: weights (samples * 1 * trials)
%
% Weights, if provided, control the size of each trial.

if nargin<2; w=[]; end

[nsamples,nchans,ntrials]=size(x);
y={};
if isempty(w)
    for k=1:ntrials
        y{k}=x(:,:,k)';
    end
else
    for k=1:ntrials
        y{k}=x(1:max(find(w(:,1,k))),:,k)';
    end
end
