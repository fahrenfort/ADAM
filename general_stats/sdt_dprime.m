function [d,beta,criterion_c, ideal_c, ideal_obs] = sdt_dprime(pHit,pFA, cormethod, nSignal, nNoise)
%  [d,beta,criterion_c, ideal_c, ideal_obs] = sdt_dprime(pHit,pFA,cormethod,nSignal,nNoise)

%  pHit         hit rate: nHits/nSignal
%  pFA          false alarm rate: nFalseAlarms/nNoise
%  cormethod    method to correct for 0 and/or 1 values
%               'arbitrary' (default) simply replacing 0 with .00001 and 1 with .99999
%               'hautus' loglinear approach corrects for extreme values (0 or 1) by adding 0.5
%               to both the number of hits and the number of false alarms of all cells, and add 1 to
%               both the number of signal trials and the number of noise trials of each cell.
%               'macmillan' Adjust only the extreme values by replacing rates of 0 with 0.5/n and
%               rates of 1 with (n-0.5)/n, where n is the number of signal or noise trials.
%               Both 'hautus' and 'macmillan' require nSignal and nNoise to be passed to the
%               function. 
%               'none' No correction for 0 and/or 1 values (resulting in Inf values)
%  nSignal      number of signal trials
%  nNoise       number of noise trials
%
%  PHIT and PFA are numerical arrays of the same shape.
%  PHIT is the proportion of "Hits":        P(Yes|Signal)
%  PFA is the proportion of "False Alarms": P(Yes|Noise)
%  All numbers involved must be between 0 and 1.
%  The function calculates the d-prime measure for each <H,FA> pair.
%  The bias/criterion values beta or c, as well as the ideal observer
%  criterion and performance scores can also be requested.
%  Requires MATLAB's Statistical Toolbox.
%
%  References:
%  * Green, D. M. & Swets, J. A. (1974). Signal Detection Theory and
%    Psychophysics (2nd Ed.). Huntington, NY: Robert Krieger Publ.Co.
%  * Macmillan, Neil A. & Creelman, C. Douglas (2005). Detection Theory:
%    A User's Guide (2nd Ed.). Lawrence Erlbaum Associates.
%  * Hautus, M. J. (1995). Corrections for extreme proportions and their biasing effects on
%    estimated values ofd′. Behavior Research Methods, Instruments, & Computers, 27(1), 46–51.
%  * Macmillan, N. A., & Kaplan, H. L. (1985). Detection theory analysis of group data: estimating
%    sensitivity from average hit and false-alarm rates. Psychological Bulletin, 98(1), 185–199.
%  * Wickens, T. D. (2002). Elementary signal detection theory (Oxford; New York). Oxford University Press.
%  
%  See also NORMINV, NORMPDF.

if nargin < 3 || isempty(cormethod)
    cormethod = 'arbitrary';
end
if nargin < 4
    nSignal = [];
end
if nargin < 5
    nNoise = [];
end
if nargin > 2
    if isempty(nSignal) || isempty(nNoise)
        error('Need to also pass the number of signal (nSignal) and noise (nNoise) trials in the dataset to be able to apply Hautus or MacMillan correction.');
    end
    nHits   = pHit .* nSignal;
    nFA     = pFA .* nNoise;
end
if strcmpi(cormethod,'arbitrary')
    pHit(pHit==1)   = .99999;
    pFA(pFA==1)     = .99999;
    pHit(pHit==0)   = .00001;
    pFA(pFA==0)     = .00001;
elseif strcmpi(cormethod,'hautus')
    nHits   = nHits + 0.5;
    nFA     = nFA + 0.5;
    pHit    = nHits./(nSignal+1);
    pFA     = nFA./(nNoise+1);
elseif strcmpi(cormethod,'macmillan')
    % replace 0 with 0.5/n and rates of 1 with (n-0.5)/n
    pHit(pHit==1)   = (nSignal(pHit==1)-0.5)./nSignal(pHit==1);
    pFA(pFA==1)     = (nNoise(pFA==1)-0.5)./nNoise(pFA==1);
    pHit(pHit==0)   = 0.5./nSignal(pHit==0);
    pFA(pFA==0)     = 0.5./nNoise(pFA==0);
end

%-- Convert to Z scores, no error checking
zHit = norminv(pHit) ;
zFA  = norminv(pFA) ;

%-- Calculate d-prime
d = zHit - zFA ;

%-- If requested, calculate BETA
if (nargout > 1)
  yHit = normpdf(zHit) ;
  yFA  = normpdf(zFA) ;
  beta = yHit ./ yFA ;
end

%-- If requested, calculate criterion c
if (nargout > 2)
  criterion_c = -(zHit + zFA)/2; % log beta, divide this by d' if you want to have a measure that is scaled by detectability
end

%-- If requested, calculate the ideal criterion and ideal observer score
%-- Pc = propS.*pHit + propN.*(1-pFA); % the overall probability of being correct
%-- Here we compute the ideal observer score, which is this probability under the ideal criterion given the ratio of signal to no signal trials
if (nargout > 3)
  propS = nSignal./(nSignal+nNoise);
  propN = 1-propS;
  ideal_c = log(propN./propS); % ideal criterion, see Wickens
  ideal_l = d./2 + ideal_c./d;
  ideal_obs = propS*(1-normcdf(ideal_l-d))+propN*normcdf(ideal_l); % ideal observer score, see Wickens
end
