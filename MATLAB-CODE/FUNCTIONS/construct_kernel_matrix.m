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

for iangle=1:length(Angles);
Angle=Angles(iangle);
    
for idep=1:nd;
    fprintf('  %.2f complete at %s\n',(idep-1)/nd*100.0,datestr(now));
    
    Scat_Depth=Scat_Depths(idep);
    
    %disp(Scat_Depth);
    
    %scale=40.0;
    basedir=sprintf('../../../../OUTPUT_FILES_%d-%1.2f-%d-%2d-DBLPERIOD/OUTPUT_FILES/',Scat_Depth*1000,Scat_Strength,Scat_Radius,Angle);

    for ista = 1:1:ns;
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
        
        %debugging
        a=max(abs(interp1(T,Daughter,KTimes)));
        b=Pmax;
        fprintf('%8e %8e %8e\n',a, b, a/b )
        
        %clf;
        %plot(P); hold on;
        %plot(S);
        
        
        
    end

end

end

save(sprintf('Kernel_Angles_x2_%1.2f_%d.mat',Scat_Strength,Scat_Radius))

end

