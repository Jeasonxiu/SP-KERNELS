%
x=gausswin(500,5);
t=1:length(x);
dt=1;
X=fftshift(fft(x));
fnyq=0.5/dt;
f=linspace(-fnyq,fnyq,length(X));

X2=X.*(2*pi*1i*f');
X3=X.*(2*pi*1i*f').^(1/2);

x2=real(ifft(ifftshift(X2)));
x3=real(ifft(ifftshift(X3)));

figure(1); clf;
hold on;
l1=plot(t,x/max(abs(x)),'-k');
l2=plot(t,x2/max(abs(x2)),'-r');
l3=plot(t,x3/max(abs(x3)),'--r');

legend([l1,l2,l3],'orig. pulse','1st deriv.','half deriv.')

xlabel('Time (s)')
ylabel('Normalized Amplitude')