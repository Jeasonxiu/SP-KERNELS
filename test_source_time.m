clear; clf;
load OUTPUT/15-Aug-2017/TEST/STFun
load OUTPUT/15-Aug-2017/TEST/I

[d]=extractSeis(Inversion.DataVector,1);

tsh=8.8;
plot(STFun.time,STFun.amplitude); hold on
plot(Inversion.kernel.KTimes+tsh,d*7.25)