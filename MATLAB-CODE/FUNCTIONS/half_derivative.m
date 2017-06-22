%
function y = half_derivative(x,dt)
X=fftshift(fft(x));
fnyq=0.5/dt;
f=linspace(-fnyq,fnyq,length(X));
X3=X.*(2*pi*1i*f').^(1/2);
y=real(ifft(ifftshift(X3)));
end