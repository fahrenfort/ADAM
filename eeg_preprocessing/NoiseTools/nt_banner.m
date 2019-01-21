function h=nt_banner(text)
%h=nt_banner(text,varargin) - annotate with text at head of figure
% 
%  h: handle to annotation
%  
%  text: string to print
% 
%  default horizontal position is 'center'
% 
% NoiseTools

c=get(gcf,'children');
y=0;

% find the top of uppermost subplot
for k=1:numel(c)
    pos=get(c(k),'position');
    y=max(y,pos(2)+pos(4));
end

 h=annotation('textbox',[0 y 1 1-y], 'linestyle', 'none', 'string', ...
     text, 'interpreter', 'none', 'verticalalignment', 'top', 'horizontalalignment', 'center');
 
 if nargout==0; clear h; end
 