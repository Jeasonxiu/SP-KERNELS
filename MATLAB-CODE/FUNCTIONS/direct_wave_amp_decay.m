function direct_wave_amp_decay(Scat_Strength,Scat_Radius)

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

Angles=15;%,30,32];

%lf=1.0/150.0;
%hf=1.0/12.0;

figure(1); clf;
hold on;

for iangle=1:length(Angles);
Angle=Angles(iangle);
    
for idep=1:nd;
    fprintf('  %.2f complete at %s\n',(idep-1)/nd*100.0,datestr(now));
    
    Scat_Depth=Scat_Depths(idep);
    
    %disp(Scat_Depth);
    
    %scale=40.0;
    basedir=sprintf('/gpfs_home/nmancine/data2/nmancine/PROJECTS/SP_RECEIVER_FUNCTIONS/KERNEL/OUTPUT_FILES_%d-%1.2f-%d-%2d/OUTPUT_FILES/',Scat_Depth*1000,Scat_Strength,Scat_Radius,Angle);

    for ista = 1:1:ns;
        station=Stations(ista);

        %load data in fast
        tmpfile=[basedir sprintf('AA.S0%03d.BXZ.semd',station)];
        %fid=fopen(tmpfile);
        %try
        %    nfoP=textscan(fid,'%f %f');
        %catch
        %    error(['No file: ' tmpfile])
        %end
        %fclose(fid);
        fid=fopen([basedir sprintf('AA.S0%03d.BXX.semd',station)]);
        nfoS=textscan(fid,'%f %f');
        fclose(fid);

        tS=nfoS{1};
        S=nfoS{2};

        Parent=S';
        

        dt=tS(2)-tS(1);
        
        [~,tmax]=max(abs(Parent));
        
  
        tmax_plus = tmax + round(0./dt);      
        maxamp=max(Parent(tmax_plus));
        
        fprintf(' %d %f\n',ista,maxamp)
        plot(ista,maxamp)
        
        
        


    end

end

end

end

