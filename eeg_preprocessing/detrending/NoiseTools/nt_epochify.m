function y=nt_epochify(x,tidx,bounds);
%y=nt_epochify(x,tidx,bounds) - extract epochs based on trigger indices
%
%  y: 3D array of epoched data (time*channels*trials)
%
%  x: raw data (time*channels)
%  tidx: array of trigger indices
%  bounds: (samples) start and stop of epoch
%  
%  if elements or tidx are fractionnary, the epochs are extracted by
%  interpolation


if nargin<3; error('!'); end

if numel(bounds)==1; bounds=[0 bounds]; end

if max(tidx)+bounds(2)>size(x,1);
    [max(tidx) bounds(2) size(x,1)]
    error('epoch stops after end of data');
end
if min(tidx)+bounds(1)<1;
    error('epoch starts before beginning of data');
end

y=zeros(bounds(2)-bounds(1)+1, size(x,2), numel(tidx));

if tidx == round(tidx)
    % integer positions
    for k=1:numel(tidx);
        y(:,:,k)=x(tidx(k)+(bounds(1):bounds(2)),:);
    end
else
    % fractionnary positions
    warning('nointeger tidx, using interpolation'); 
    for k=1:numel(tidx);
        y(:,:,k)=interp1(1:size(x,1),x,tidx(k)+(bounds(1):bounds(2)));
    end
end