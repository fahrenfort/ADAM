function nt_linestyles(h,property,values)
%nt_stylelines(h,property,values) - apply different styles to lines of plot
%
%  h: handle to plot (default:gca)
%  property: property to set (default:linewidth)
%  values: cell array of values, one for each line
%
% Values may be a numerical matrix, in which case each row is used as a
% value.
%
% Styles are applied to children of h in reverse order (ie in order of plot
% commands).  May produce unexpected results if there are childern other
% than plot lines.
% 
% NoiseTools
%
%{ 
examples of properties:
                 Color: [0 0 1]
             LineStyle: '-'
             LineWidth: 0.5000
                Marker: 'none'
            MarkerSize: 6
       MarkerEdgeColor: 'auto'
       MarkerFaceColor: 'none'
%}

if nargin<1 || isempty(h); h=gca; end
if nargin<2; property=[]; end
if nargin<3; error('!'); end

if isempty(property); property='linewidth'; end

c=get(h,'children');

if ~isa(values, 'cell'); 
    values=num2cell(values,1);
end


for k=1:numel(c);
    set(c(numel(c)-k+1),property,values{mod(k-1,numel(c))+1})
end

