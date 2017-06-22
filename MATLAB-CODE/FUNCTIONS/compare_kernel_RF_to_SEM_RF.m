function compare_kernel_RF_to_SEM_RF()
    
    [~,seis1]=     Kernel_Summation(1,2,1,0.02,2500,0,150,160);
    [~,seis2]=     Kernel_Summation(1,2,1,0.02,2500,0,150,160);
    [KTimes,seis3]=Kernel_Summation(1,2,1,0.02,2500,0,150,160);
   
    nf=max(abs(seis1(KTimes>-10)));
    l1=plot(KTimes,seis1/nf,'--k','LineWidth',2);
    hold on
    nf=max(abs(seis2(KTimes>-10)));
    l2=plot(KTimes,seis2/nf,'--g');
    nf=max(abs(seis3(KTimes>-10)));
    l3=plot(KTimes,seis3/nf,'--m');
    ylabel('RF Amplitude')
    xlabel('Time (s)')
    
%     for iseis = [50:70];
%         
%         Daughter=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/CASE4.BAK.OCT_6/OUTPUT_FILES_50000-40000-1/OUTPUT_FILES/AA.S0%03d.BXP.semd',iseis));
%         Parent=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/CASE4.BAK.OCT_6/OUTPUT_FILES_50000-40000-1/OUTPUT_FILES/AA.S0%03d.BXS.semd',iseis));
%         
%         %dt=Parent(2,1)-Parent(1,1);
%         %Parent(:,2)=bpfilt(Parent(:,2)',dt,0.01,0.05);
%         %Daughter(:,2)=bpfilt(Daughter(:,2)',dt,0.01,0.05);
%         
%         [~,imax]=max(abs(Parent(:,2)));
%         time=Parent(:,1) - Parent(imax,1);
%         dtmp=interp1(time,Daughter(:,2),KTimes)';
% 
%         nf=max(abs(dtmp));
%         l4=plot(KTimes,dtmp/nf,'-r');
%     end
    
%     for iseis = [50:70];
%         Daughter=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/CASE4.BAK.OCT_5/OUTPUT_FILES_50000-40000-1/OUTPUT_FILES/AA.S0%03d.BXP.semd',iseis));
%         Parent=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/CASE4.BAK.OCT_5/OUTPUT_FILES_50000-40000-1/OUTPUT_FILES/AA.S0%03d.BXS.semd',iseis));
%         [~,imax]=max(abs(Parent(:,2)));
%         time=Parent(:,1) - Parent(imax,1);
%         dtmp=interp1(time,Daughter(:,2),KTimes)';
% 
%         nf=max(abs(dtmp));
%         l4=plot(KTimes,dtmp/nf,'-c');
%     end
    
    for iseis = [50:70];
        Daughter=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/KERNEL_GMSH/KERNEL/TWOLAY.4s/OUTPUT_FILES_0.05-23/OUTPUT_FILES/AA.S0%03d.BXP.semd',iseis));
        Parent=load(sprintf('~/data2/nmancine/nmancine/SP-RF-SYN/KERNEL_GMSH/KERNEL/TWOLAY.4s/OUTPUT_FILES_0.05-23/OUTPUT_FILES/AA.S0%03d.BXS.semd',iseis));
        [~,imax]=max(abs(Parent(:,2)));
        time=Parent(:,1) - Parent(imax,1);
        dtmp=interp1(time,Daughter(:,2),KTimes)';

        nf=max(abs(dtmp));
        l5=plot(KTimes,dtmp/nf,'-m');
    end

    legend([l1,l5],'Kernel summation 1','SEM calc. + mesh','Location','SouthWest')
    uistack(l1,'top')
    
end
    
 function [KTimes,seis]=Kernel_Summation(tDeci,xDeci,zDeci,Scat_Strength,Scat_Radius,stamin,stamax,interface_depth)
    tmp=sprintf('Kernel_Angles_x2_%1.2f_%d.mat',Scat_Strength,Scat_Radius);
    load(tmp,'Kernel','Stations','Scat_Depths','KTimes');
  
    ip=1; %just use first angle calc for now
    
    Kernel=squeeze(Kernel(1:zDeci:end,1:xDeci:end,ip,1:tDeci:end));
    Stations=Stations(1:xDeci:end);
    Scat_Depths=Scat_Depths(1:zDeci:end);
    KTimes=KTimes(1:tDeci:end);
    nTimes=length(KTimes);
    
    seis=zeros(length(KTimes),1);
    figure(1);clf;
    
    for iscat =1:length(Scat_Depths);
        for ioff =1:length(Stations);
            tmp=squeeze(Kernel(iscat,ioff,:));
            if Stations(ioff) < stamin || Stations(ioff)>stamax; 
                continue
            end
            if Scat_Depths(iscat)<=interface_depth;
                seis=seis-1.0*tmp;
            else
                seis=seis+0.0*tmp;
            end
        end
    end
    
    tmp=diff(Scat_Depths);
    fprintf('dz = %f\n',tmp(1))
    
 end
 
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
