function [y] = bpfilt(x,dt,lf,hf)

% ********* Function Description *********
%
% Bandpass filter a time seriers.
%
% [Y] = bpfilt(X,DT,LF,HF)
%
% Take a time series, X, sampled at DT and
% filter it with a 2nd order, 2 pass
% butterworth filter between frequencies
% LF and HF. If X is a matrix, this will
% filter the individual rows of X.
%
% ****************************************
% *                                      *
% *  Written by David L. Abt - May 2008  *
% *                                      *
% *  Taken from code written by Michael  *
% *  Bostock and Stephane Rondenay, and  *
% *  used by Kate Rychert.               *
% *                                      *
% *  Email: David_Abt@brown.edu          *
% *                                      *
% ****************************************

nyq     = 0.5/dt;           % Nyquist Frequency
wn      = [lf/nyq,hf/nyq];
[b,a]   = butter(2,wn);
for ix=1:length(x(:,1))
    y(ix,:) = filtfilt(b,a,double(x(ix,:)));  % Edited by Ved because waveforms
    % are single not double, but filtfilt
    % requires double
end

end