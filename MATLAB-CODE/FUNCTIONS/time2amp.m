function [A] = time2amp(DT,trial_times,amps)

%Simple gaussian 
extrapval=0.0;

A=interp1(trial_times,amps,DT,'pchip',extrapval);

%A = exp(-DT.^2./(2.0*tchar^2));

%Ricker wavelet
%A = (1-(DT./tchar).^2) .* A;

%Deriv of Ricker wavelet
%A = -2*(tchar-DT)/tchar^2 .* A + (1-(DT/tchar).^2) .* (-DT/tchar^2) .* A;

%A = -A .* (DT.^3 - 3*tchar^2 *DT)/tchar^4;

end