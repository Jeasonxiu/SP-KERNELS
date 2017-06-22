function [out]=fractional_deriv(a,n,t)

fun = @(w) 1/(2*pi)^0.5 * a*exp(-1/2 * a^2 * w.^2) .* exp(1i*w*t) .*(1i*w).^n;
amp=integral(@(w)fun(w),-Inf,Inf);

out=real(amp);

end


