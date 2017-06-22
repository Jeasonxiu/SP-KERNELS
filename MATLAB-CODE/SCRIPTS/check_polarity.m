clear all;
figure(1)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 800, 1000]); %<- Set size

Scat_Depth=100.0;
scale=2000.0;
addpath('~/PROJECTS/ARRAY_STACK/ReceiverFunctions/Make_Receiver_Functions/CLEAN/Scattered_Waves/ETMTM/')
s1=1;
s2=5;
s3=75;
Stations=s1:s2:s3;
ns=length(Stations);

lf=1.0/150.0;
hf=1.0/12.0;

for ipol=1:2;
    if ipol==1;
        pol='';
        linecolor='b';
    else
        pol='-';
        linecolor='r';
    end
    
    basedir=sprintf('/Users/mancinelli/PROJECTS/SP-RF-SYN/KERNEL/OUTPUT_FILES_%d000-%s0.01-1/OUTPUT_FILES/',Scat_Depth,pol);


for ista = 1:ns;
    station=Stations(ista);

    %load data in fast
    fid=fopen([basedir sprintf('AA.S00%02d.BXP.semd',station)]);
    nfoP=textscan(fid,'%f %f');
    fclose(fid);
    fid=fopen([basedir sprintf('AA.S00%02d.BXS.semd',station)]);
    nfoS=textscan(fid,'%f %f');
    fclose(fid);
    tP=nfoP{1};
    P=nfoP{2};
    tS=nfoS{1};
    S=nfoS{2};

    Parent=S';
    Daughter=P';
    TB=4;
    NT=3;
    dt=tS(2)-tS(1);
    
    Parent=bpfilt(Parent,dt,lf,hf);
    Daughter=bpfilt(Daughter,dt,lf,hf);

    %plot(tS,P);
    %pause;

    [~,tmax]=max(abs(Parent));

    tmax_plus = tmax + round(25./dt);
    tmax_minus = tmax - round(50./dt);

    if tmax_minus<1;
        tmax_minus=1;
    end
    
    fprintf('%d %d %f s\n',tmax_plus,length(Parent),(length(Parent)-tmax_plus)*dt)

    %Center parent phase at 25s from the window end
    Parent=Parent(tmax_minus:tmax_plus);
    Daughter=Daughter(tmax_minus:tmax_plus);
    T=tS(tmax_minus:tmax_plus);

    %Add noise 
    %Daughter=Daughter + randn(1,length(Daughter));

    %For ETMTM
    %[TRF,RF] = ETMTM(Parent,Daughter,TB,NT,'synth',dt);
    
    %For ITDD
    Time=(-50:dt:5);
    nIter=10;
    PulseFreq=0.5;
    [RF] = ITDD(Parent,Daughter,dt,Time,'Sp',nIter,PulseFreq);
    TRF=Time;

    % Trim RF
    RF=RF(TRF>-50);
    TRF=TRF(TRF>-50);

    %plot(TRF,RF*scale+station,linecolor)
    subplot(1,3,1)
    plot(T,Daughter*scale+station,linecolor)
    hold on;
    subplot(1,3,2)
    plot(T,Parent*5.0+station,linecolor)
    hold on;

    subplot(1,3,3)
    plot(TRF,RF*scale*30.0+station,linecolor)
    hold on;

end
end

subplot(1,3,1)
title(sprintf('Scatterer at a Depth of %d km',Scat_Depth))
ylabel('Parent Waveforms')
xlabel('Time')

subplot(1,3,2)
title(sprintf('Scatterer at a Depth of %d km',Scat_Depth))
ylabel('Daughter Waveforms')
xlabel('Time')

subplot(1,3,3)
title(sprintf('Scatterer at a Depth of %d km',Scat_Depth))
ylabel('RF')
xlabel('Time')

%pause;
%print('polarity','-dpng','-r300');
print('polarity','-depsc2','-r300');

