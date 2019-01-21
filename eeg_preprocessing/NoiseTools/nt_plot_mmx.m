function nt_plot_mmx(a,b,c)
%nt_plot_mmx - plot data using min-max pairs 
%
%   nt_plot_mmx(x): plot x
%   nt_plot_mmx(abscissa,x): plot x using abscissa
%   nt_plot_mmx(abscissa,x,N): plot x using abscissa, with specified target number of min-max pairs [default: 1000]

error('not yet implemented');

if nargin < 3 || isempty(c); N=1000; else N=c; end
if nargin < 2 || isempty(b); 
    x=a;
    abscissa=1:size(x,1); 
else
    x=b;
    abscissa=a; 
end

[mmx_pairs,a]=nt_mmx(x);
abscissa=abscissa(a+1);
plot(abscissa,x);


