% Testing isochrons
%   The P-to-S isochron is a parabola
%   The S-to-P isochron is a hyperbola

t0=10;
vp=4.5;
vs=8.0;
angle=24;

x0s=-600:600;

% ys=sqrt( (t0*vp + -x0s * sind(angle) *vp/vs).^2 + x0s.^2   );
% 
% plot(x0s,-ys)

A=vp*t0;
B=vp/vs;

NUM=x0s.^2-A.^2;
DEN=B.^2 - 1;

ys=sqrt( NUM / DEN);

plot(x0s,-ys)
