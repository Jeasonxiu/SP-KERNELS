function [T] = timeshifts(Tdirect,Tscat,xs,zs)
T10=Tdirect;
T20=Tscat;

[XS0, ZS0] = meshgrid(xs,zs);

tref=interp2(XS0,ZS0,T10,0,0);

T = T10'+T20 - tref;
end