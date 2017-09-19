%
%
clear;
figure(1)
alpha=7.920;
beta=4.400;
rho=3.300;
dbeta_over_beta=1;

%from scat pattern
fac = rho;
fac = fac/sqrt(alpha); %my empirical correction

%from ray tracing (eqivalent to fac * sqrt(1/2*pi*r))
load geom_halfspace
subplot(2,2,1)
A=amp'*fac;
pcolor(xs,zs,log10(A)); shading interp; colorbar;
title('Ray Tracing (A)')

%from Junlin
[XS,ZS]=meshgrid(xs,zs);
R=sqrt(XS.^2 + ZS.^2);
G=1/(4*alpha^2)*sqrt(2*alpha./(pi*R));
B=G;
subplot(2,2,2)
pcolor(xs,zs,log10(B)); shading interp; colorbar;
title('Junlin Greens (B)')

%from Shragge
G=(1/4*alpha)*sqrt(2/pi/alpha/alpha);
G=G./R;
subplot(2,2,3)
pcolor(xs,zs,log10(G)); shading interp; colorbar;
title('Shragge')

%from Hudson
[XS,ZS]=meshgrid(xs,zs);
R=sqrt(XS.^2 + ZS.^2);
G=rho*sqrt(1/2/pi/alpha./R);
B=G;
subplot(2,2,4)
pcolor(xs,zs,log10(B)); shading interp; colorbar;
title('Hudson (B)')


figure(2)
subplot(3,1,1)
TMP=abs(dpdx)'./cosd(THE2).^2*alpha;
A=sqrt(TMP');
pcolor(xs,zs,A);shading interp; colorbar;
title('Ray Tracing')
subplot(3,1,2);
B=1./sqrt(R);
%B=R;
pcolor(xs,zs,B); shading interp; colorbar;
title('$1/\sqrt{r}$','interpreter','latex')
subplot(3,1,3)
contourf(xs,zs,log10(A./B)); shading interp; colorbar;
caxis([-.05,.05])
disp(nanmedian(nanmedian(A./B)));