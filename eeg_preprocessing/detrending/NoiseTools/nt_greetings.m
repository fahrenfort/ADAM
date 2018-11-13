function nt_greetings(reset)
%nt_greetings - display message the first time the toolbox is used

persistent nt_greeted
if nargin>0; nt_greeted=reset; return; end

if isempty(nt_greeted)
    
    display(' ')
    display(['NoiseTools, version ',nt_version]);
    display('http://audition.ens.fr/adc/NoiseTools');
    display('Please cite relevant methods papers.');
    display(' ');
    
    nt_greeted=1;
    
end