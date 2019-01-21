function nt_colorlines(h,permutation)
%nt_colorlines(h,permutation) - apply different colors to lines of plot
%
%  h: handle to plot (default:gca)
%  permute: permutation to apply to colors
%
% Colors are applied to children of h in reverse order (ie in order of plot
% commands).  May produce unexpected results if there are childern other
% than plot lines.\
% 
% NoiseTools
% See nt_stylelines, nt_widthlines.

if nargin<1 || isempty(h); h=gca; end
if nargin<2; permutation=[]; end

colororder=get(h,'colororder');
if ~isempty(permutation); 
    colororder=colororder(permutation,:);
end
c=get(h,'children');

for k=1:numel(c);
    set(c(numel(c)-k+1),'color',colororder(1+mod(k-1,size(colororder,1)),:))
end

