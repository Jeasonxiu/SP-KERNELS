classdef GaussianFilter
   properties
   end
   methods (Static)
       function tr_fil = filter_trace(tr,stdev,dt)
           npts=round(stdev/dt*3.0);
           alpha=1/stdev*dt*npts;
           %fprintf('npts, alpha = %d, %f\n',npts,alpha)
           gw=gausswin(npts,alpha);
           %wvtool(gw)
           tr_fil=conv(tr,gw,'same')/sum(gw);
       end
 
       function test
           function [tr, time] = generate_test_trace
                dt=1.0;
                time = -500:dt:500;
                tr=gausswin(length(time),30);
           end
           [tr,time]=generate_test_trace;
           stdev=50;
           trfil=GaussianFilter.filter_trace(tr,stdev,dt);
           plot(time,tr); hold on
           plot(time,trfil)
       end
   end
   
end