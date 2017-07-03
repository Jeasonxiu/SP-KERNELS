function [amp] = geom_spreading(Pscat,xs,zs,model)
%%
P = Pscat / 111.1;

%[vps,~]=get_velocity_from_profile('MIGRA/myvmod.nd',zs);
[vps,~]=get_v(model,zs);

[VPS,~]=meshgrid(vps,xs);

U1=1./VPS;
THE1=asind(P./U1);
THE2=asind(P/U1(1,1));

%fprintf('surface velocity = %f\n',1/U1(1,1));

%take derivative and resample
tmp=diff(P);
xmid=(xs(1:(end-1))+xs(2:end))/2;
dpdx=interp2(xmid,zs,tmp',xs,zs);


%This is modified from Eq 6.23 in Shearers Intro to
%Seismology (derived for 2-D)
energy=(1/2/pi)*abs(dpdx)'./cosd(THE1)./cosd(THE2)./U1;

amp=(sqrt(energy));

%figure(1)
%subplot(2,1,1)
%
%contourf(xs,-zs,log10(amp')); shading flat; colorbar
%
%subplot(2,1,2)
%
%[XS,ZS]=meshgrid(xs,zs);
%D=sqrt(XS.^2 + ZS.^2);
%contourf(xs,-zs,log10(1./sqrt(D))); shading flat; colorbar
%%
end