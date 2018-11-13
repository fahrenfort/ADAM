function w=nt_growmask(w,margin)
%ww=nt_growmask(w,margin) - widen mask
%
%  ww: widened mask (samples * 1 or samples * 1 * trials)
%
%  w: mask
%  margin: samples, number of samples by which to widen the mask
%
% NoiseTools

[m,n,o]=size(w);
w=(w==0);
w=nt_unfold(w);
w=filter(ones(margin+1,1),1,w)>1;
w=flipud(filter(ones(margin+1,1),1,flipud(w))>1);
w=(w==0);
w=nt_fold(w,m);