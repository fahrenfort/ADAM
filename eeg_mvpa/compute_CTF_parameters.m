function avgctfstruct = compute_CTF_parameters(avgctfstruct,gsettings)
% function new_avgctfstruct = compute_CTF_parameters(avgctfstruct,gsettings)
% computes CTF parameters sigma (FWHM), amplitude, as well as linear slope
% for for all subjects and for the group average. The sigma from the group
% average can be input as parameter into the forward encoding model
% functions to achieve an iterative approach (i.e. fitting the initial
% model using a delta function and then using the sigma from this function
% as input for a new forward encoding model fit.
% Takes as input the output from plot_CTF.m with the '2D' option (which
% resulted in the average CTF over a set time window)
% All CTFs are first made symmetrical by averaging equidistant points from
% the horizontal midline.
% Currently, the are not normalized to have 0 as minimum and sum/AUC 1
% because I am not sure what the theoretical implications are of doing this
% If you want to suppress plotting the results set 
% gsettings.plotfittype ='', otherwise set to 'gaus' or 'slope' 

% J.J.Fahrenfort, VU 2016

% getting relevant settings
plottype = '3D';
plotfittype = 'slope'; % can be 'gaus' or 'slope'
v2struct(gsettings);

% only works on 2D
if ~strcmpi(plottype,'2D')
    disp('currently only computes CTF parameters for 2D data, easy to extend for future versions');
    return;
end

% figure for group average
nStats = numel(avgctfstruct);
if ~isempty(plotfittype); 
    fh = figure('Name','group average CTF fits','NumberTitle','off');
    set(fh,'color','w');
    if nStats > 1
        set(fh, 'Position', get(0,'Screensize'));
    end
end

% loop around stats for group average
for cStats = 1:nStats
    % get
    localstruct = avgctfstruct(cStats);
    % unpack
    v2struct(localstruct);
    % compute average
    avgCTF = mean(indivCTFs,1);
    if ~isempty(plotfittype)
        subplot(numSubplots(nStats,1),numSubplots(nStats,2),cStats);
    end
    % compute relevant parameters
    [gaus_fit, gaus_real, sigma, A] = fit_ctf(avgCTF);
    [slope_fit, slope_real, slope] = fit_slope(avgCTF);
    avgctfstruct(cStats).basis_set = gaus_fit;
    avgctfstruct(cStats).sigma = sigma;
    avgctfstruct(cStats).A = A;
    avgctfstruct(cStats).slope = slope;
    % plot
    if strcmpi(plotfittype,'gaus')
        plot_fit(gaus_fit,gaus_real,condname);
    elseif strcmpi(plotfittype,'slope')
        plot_fit(slope_fit,slope_real,condname);
    end
end

% loop around stats and loop around individual subjects (not the most efficient code but it works and easy to read)
for cStats = 1:nStats
    % get
    localstruct = avgctfstruct(cStats);
    % unpack
    v2struct(localstruct);
    % make figure
    if ~isempty(plotfittype); 
        fh = figure('Name',['individual subject CTF fits for ' condname],'NumberTitle','off');
        set(fh,'color','w');
        if nStats > 1
            set(fh, 'Position', get(0,'Screensize'));
        end
    end
    nSubj = size(indivCTFs,1);
    for cSubj = 1:nSubj
        % compute average
        subjCTF = squeeze(indivCTFs(cSubj,:));
        if ~isempty(plotfittype); 
            subplot(numSubplots(nSubj,1),numSubplots(nSubj,2),cSubj);
        end
        % compute relevant parameters
        [gaus_fit, gaus_real, sigma, A] = fit_ctf(subjCTF);
        [slope_fit, slope_real, slope] = fit_slope(subjCTF);
        avgctfstruct(cStats).indivBasis_set(cSubj,:) = gaus_real;
        avgctfstruct(cStats).indivSigma(cSubj) = sigma;
        avgctfstruct(cStats).indivA(cSubj) = A;
        avgctfstruct(cStats).indivAUC(cSubj) = sum(subjCTF-min(subjCTF)); % normalized to the minimum
        avgctfstruct(cStats).indivSlope(cSubj) = slope;
        % plot
        if strcmpi(plotfittype,'gaus')
            plot_fit(gaus_fit,gaus_real,cSubj);
        elseif strcmpi(plotfittype,'slope')
            plot_fit(slope_fit,slope_real,cSubj);
        end
    end
end



