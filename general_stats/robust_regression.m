function [rsquare, stats, beta] = robust_regression(predictor, data)
% compute R2 and beta
[beta, stats] = robustfit(predictor,data);
sse = stats.dfe * stats.robust_s^2;
phat = beta(2)*data + beta(1);
ssr = norm(phat-mean(phat))^2;
rsquare = 1 - sse / (sse + ssr);