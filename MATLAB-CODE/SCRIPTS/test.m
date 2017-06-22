clear; clf;
scale=4.0;
curdir='/Users/mancinelli/PROJECTS/SP-RF-SYN/KERNEL';
basedir='/Users/mancinelli/PROJECTS/SP-RF-SYN/KERNEL/OUTPUT_FILES_150000-0.10-1/OUTPUT_FILES/';
cd(basedir)
addpath('~/PROJECTS/ARRAY_STACK/ReceiverFunctions/Make_Receiver_Functions/CLEAN/Scattered_Waves/ETMTM/')

for station = 50;
    
    tmp=sprintf('AA.S00%2d.BXP.semd',station);
    nfoP=load(sprintf('AA.S00%02d.BXP.semd',station));
    nfoS=load(sprintf('AA.S00%02d.BXS.semd',station));

    tP=nfoP(:,1);
    P=nfoP(:,2)*scale*25.0+station;

    tS=nfoS(:,1);
    S=nfoS(:,2)*scale+station;

    Parent=S';
    Daughter=P';
    TB=4;
    NT=7;
    dt=tS(2)-tS(1);
    
    [~,tmax]=max(abs(Parent));
 
    tmax_plus = tmax + round(25./dt);
    
    Parent=Parent(1:tmax_plus);
    Daughter=Daughter(1:tmax_plus);
    T=tS(1:tmax_plus);
    
    
    %Debug;
    Daughter=Daughter*0.0;
    Parent=Parent*0.0;
    
    Daughter(tmax-1000.0)=1000.0;
    Parent(tmax)=1000.0;
    
    Daughter=Daughter + randn(1,length(Daughter));

    [TRF,RF] = ETMTM(Parent,Daughter,TB,NT,'synth',dt);
    plot(TRF,RF,'r')
    hold on;
    %plot(TRF,RF*100.0+station,'black');
    
    %plot(T,Parent+station)

end
cd(curdir)





