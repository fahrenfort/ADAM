% DPRIME  - Signal-detection theory sensitivity measure.
%
%  d = dprime(pHit,pFA)
%  [d,beta] = dprime(pHit,pFA)
%
%  PHIT and PFA are numerical arrays of the same shape.
%  PHIT is the proportion of "Hits":        P(Yes|Signal)
%  PFA is the proportion of "False Alarms": P(Yes|Noise)
%  All numbers involved must be between 0 and 1.
%  The function calculates the d-prime measure for each <H,FA> pair.
%  The criterion value BETA can also be requested.
%  Requires MATLAB's Statistical Toolbox.
%
%  References:
%  * Green, D. M. & Swets, J. A. (1974). Signal Detection Theory and
%    Psychophysics (2nd Ed.). Huntington, NY: Robert Krieger Publ.Co.
%  * Macmillan, Neil A. & Creelman, C. Douglas (2005). Detection Theory:
%    A User's Guide (2nd Ed.). Lawrence Erlbaum Associates.
%  
%  See also NORMINV, NORMPDF.

function [d,beta] = sdt_dprime(pHit,pFA)

if pHit == 1
    pHit = .99999;
end
if pHit == 0
    pHit = .00001;
end
if pFA == 0
    pFA = .00001;
end
if pFA == 1
    pFA = .99999;
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

%%  Return DPRIME and possibly BETA
%%%%%% End of file DPRIME.M
