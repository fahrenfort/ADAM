% test nt_bias_fft

nsamples=10000;
nchans=10;
SNR=0.0001;
f=20; % Hz
sr=1000;

if 0
x=randn(nsamples,nchans);
x(:,1)=SNR*sin(2*pi*f*(1:nsamples)'/sr);
x=x*randn(nchans);
end

figure(1); clf
nt_spect_plot2(x,512,[],[],sr);


[c0,c1]=nt_bias_fft(x,[20]/sr,2048); [todss,pwr0,pwr1]=nt_dss0(c0,c1); 
figure(2); clf;
plot(pwr1./pwr0,'.-')
z=nt_mmat(x,todss);
figure(3); clf
nt_spect_plot2(z,512,[],[],sr);