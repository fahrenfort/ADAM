function x=nt_dsample(x,factor)
%y=nt_dsample(x,factor) - downsample by averaging neighboring samples
%
%  y: downsampled data
% 
%  x: data to downsample (2 or 3D)
%  factor: downsampling factor (can be fractional)
%
% Downsampling is performed along columns.  If size(x,1) is not multiple of
% factor, it is truncated.
% 
% The data are lowpass filtered by convolution with a square window, which
% ensures minimal temporal distortion of the waveform. However it does not
% strongly attenuate frequency components beyond the Nyquist frequency, so
% it is not optimal from a frequency-domain point of view. If this is a
% concern, uses nt_resample() instead.
%
% NoiseTools

if nargin<2; error('!'); end
if factor==1; return; end
if factor ~= round(factor); error('factor must be integer'); end

if ndims(x)>2;
    d=size(x);
    x=reshape(x,[d(1),prod(d(2:end))]);
    x=nt_dsample(x,factor);
    x=reshape(x,[size(x,1),d(2:end)]);
    return
end

[m,n]=size(x);
a=floor(m/factor);
b=rem(m,factor);

if b>0;
    xx=x((a*factor+1):end,:); % extra bit
    x=x(1:a*factor,:);
end

x=reshape(x,[factor,a,n]);
x=mean(x,1);
x=shiftdim(x,1);

% buggy, dunno why, simpler to remove
% if b>0
%     x=[x;mean(xx,1)]; 
% end

