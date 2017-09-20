function construct_kernel_matrix(Scat_Strength,Scat_Radius)
%Constructs kernel matrix from SEM seismograms.
% N.B. Data intesive -- works best if run locally.
% N.J. Mancinelli -- June 2017
%
%
%

d1=0;
d2=5;
d3=300;
Scat_Depths=d1:d2:d3;
nd=length(Scat_Depths);

s1=1;
s2=1;
s3=150;
Stations=s1:s2:s3;
ns=length(Stations);

Angles=[15,20,25];

%start earlier to avoid S-to-S scattering
KTimes=-0:-0.2:-30;
Kernel=zeros(nd,ns,length(Angles),length(KTimes));

%lf=1.0/150.0;
%hf=1.0/12.0;

counter=0;
totcount=length(Angles)*nd*ns;
%tstart=datetime;
%tlast=datetime;

for iangle=1:length(Angles);
Angle=Angles(iangle);
    
for idep=1:nd;
    
    
    Scat_Depth=Scat_Depths(idep);
    
    %disp(Scat_Depth);
    
    %scale=40.0;
    basedir=sprintf('KERNEL-SEM/OUTPUT_FILES_%d-%1.2f-%d-%2d/OUTPUT_FILES/',Scat_Depth*1000,Scat_Strength,Scat_Radius,Angle);

    for ista = 1:1:ns;
        
        counter=counter+1;
        %telap=datetime-tstart;
        
        %tinc=datetime-tlast;
        %tlast=datetime;
        %trem=tinc*(totcount-counter);
        imod=mod(counter,20);
        %trem_vector(imod+1)=trem;
        if imod==0;
            %plot(telap,mean(trem_vector),'^r'); hold on
            fprintf('  %.5f %%.\n',counter/totcount*100);
        end
        
        station=Stations(ista);

        %load data in fast
        tmpfile=[basedir sprintf('AA.S0%03d.BXP.semd',station)];
        fid=fopen(tmpfile);
        try
            nfoP=textscan(fid,'%f %f');
        catch
            error(['No file: ' tmpfile])
        end
        fclose(fid);
        fid=fopen([basedir sprintf('AA.S0%03d.BXS.semd',station)]);
        nfoS=textscan(fid,'%f %f');
        fclose(fid);
        tP=nfoP{1};
        P=nfoP{2};
        tS=nfoS{1};
        S=nfoS{2};

        Parent=S';
        Daughter=P';

        dt=tS(2)-tS(1);
        
        [Pmax,tmax]=max(abs(Parent));
        tdif=tP(tmax);
  
        tmax_plus = tmax + round(0./dt);
        tmax_minus = tmax - round(50./dt);
        
        if tmax_minus<1;
            tmax_minus=1;
        end
        

        
        %%
        %figure(1)
        %       
        %clf; plot(tP,P); hold on; plot(tS,S,'-r')
        %xlim([tP(tmax_minus),tP(tmax_plus)])
        %
        %
        %%
        
        Daughter=Daughter(tmax_minus:tmax_plus);
        T=tS(tmax_minus:tmax_plus)-tdif;
        Kernel(idep,ista,iangle,:)=interp1(T,Daughter,KTimes);
        
        
        %clf;
        %plot(P); hold on;
        %plot(S);
        
        
        
    end

end

end

save(sprintf('Kernel_Angles_x2_%1.2f_%d_new.mat',Scat_Strength,Scat_Radius))

end

