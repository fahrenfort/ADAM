function plot_fit(fit_data,real_data,condname, gofcriterion)
% plots real data, fitted data and evaluates fit using goodness of fit
% criterion
% J.J.Fahrenfort, VU 2016

if nargin < 4
    gofcriterion = .5;
end

% plot actual
y = real_data;
x = 1:numel(y);
plot(x,y);

% good fit?
GOF = goodnessOfFit(fit_data',real_data','NRMSE');

% plot fitted data
hold on;
y = fit_data;
plot(x,y,'r');
axis square;
set(gca,'FontSize',22);
sameaxes('xyzc',gcf());
if isnumeric(condname)
    ntitle(sprintf('subj%s, GOF = %0.2f',num2str(condname),GOF));
    if GOF < gofcriterion
        text(1,.1,'bad fit, throw out?','Color','red');
    end
else
    title(sprintf('%s, GOF = %0.2f',regexprep(condname,'_',' '),GOF));
    legend({'real data','fitted data'});
end