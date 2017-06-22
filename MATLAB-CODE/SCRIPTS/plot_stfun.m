a=0.8;
n=3.0;

times=-10:0.1:10;
amps=zeros(length(times),1);

for itime=1:length(times)
   time=times(itime);
   amps(itime)=fractional_deriv(a,n,time);
end

figure(1)
plot(times,amps)
    