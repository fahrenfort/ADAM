function res = jackT(x,y)
%res = jackT(x) or res = jackT(x,y)
%x and y are vectors containg n leave-one-out estimates, plus the 
%grand-average estimate in the last row
%
%written by Heinrich René Liesefeld, September 2017
%For details and formula, see:
%Miller, J.O., Patterson, T.,& Ulrich, R. (1998). Jackknife-based method
%    for measuring LRP onset latency differences. Psychophysiology, 35,
%    99-115. doi:10.1111/1469-8986.3510099

if exist('y', 'var')
    x = x - y;
end
n = length(x) - 1;
res.mean  = mean(x(1:end-1));
res.sd    = sqrt(((n - 1) / n) * sum((x(1:end-1) - res.mean) .^2));
res.tstat = x(end) / res.sd; %use the estimate from the full GA instead of the mean of the leave-one-outs
res.df    = n  -1;
res.p     = 2 * tcdf(-abs(res.tstat), res.df);
end