function back_project_synthetics(InversionParams)
% 
% Inputs: Inversion Params
%
% Files written: Mat file containing model and relevant metadata
%

    skipSta=InversionParams.skipSta;
    tDeci=InversionParams.tDeci;
    dxin=InversionParams.dxin;
    dzin=InversionParams.dzin;
    Nus=InversionParams.Nus;
    Norm_Opts=InversionParams.Norm_Opts;
    nIterMax=InversionParams.nIterMax;
    direction=InversionParams.direction;
    Kernel_Type=InversionParams.Kernel_Type;
    ImagingMethod=InversionParams.ImagingMethod;
    saveFilename=InversionParams.saveFilename;

    if (ImagingMethod == 1)
        fprintf('Performing back-projection migration\n')
        %Thin kernels
        tchar=0.1;
        nderiv=0.0;
        %One iteration
        nIterMax=1;
        %Don't take data difference
        TakeDifferences=false;
        %Don't regularize
        Norm_Opts=0;
        %Deconvolve
        DeconvolveParentWaveform=true;
        
    elseif (ImagingMethod == 2)
        fprintf('Performing C-G Inversion\n')
        %FF kernels
        tchar=0.8;
        nderiv=3.0;
        %Many iterations
        %use nIterMax
        %Take data difference
        TakeDifferences=true;
        %Don't regularize
        Norm_Opts=2;
        %Don't deconvolve
        DeconvolveParentWaveform=false;
    else
        return
    end

    if Kernel_Type == 1;
        path='/Volumes/nmancine/data2/nmancine/PROJECTS/SP_RECEIVER_FUNCTIONS/KERNEL/SP-KERNELS/DATA';
        [Kernel,X,Y,Z,KTimes,nTimes] = load_kernel(path,tDeci);
    elseif Kernel_Type ==2;
        [Kernel,zs,xs,~,KTimes] = analytical_kernel_halfspace(26);
        nTimes=length(KTimes);
        [X,Y,Z]=meshgrid(xs,zs,1:nTimes);
        TMP(:,:,:)=Kernel(:,:,1,:);
        Kernel=TMP;
    elseif Kernel_Type ==3;
        [~,betatmp]=get_velocity_from_profile('MIGRA/myvmod.nd',300);
        incangs=[20,23,26];
        rps=sind(incangs)/betatmp;
        %fprintf('Debug 1: %f %f %f\n',incangs(1),betatmp,rps)
        fprintf('Computing analytical kernel... ')
        %[Kernel,X,Y,Z,KTimes,nTimes] = analytical_kernel_layered(rps,tchar,nderiv);
        model.hs=[60.0,60.0,180.0];
        model.vp=[5.660,7.920,7.280];
        model.vs=[3.2,4.4,4.05];
        xs=(-600:1:150)';
        zs=0:1:300;
        [Kernel,X,Y,Z,KTimes,nTimes] = analytical_kernel_layered(rps, xs, zs, tchar, nderiv, model);
        fprintf('done.\n')
    else
       error('Bad Kernel_Type: %d',Kernel_Type);
    end
    
    Kernel=real(Kernel);
    Kernel(isnan(Kernel))=0;
    
    %Set model parameters
    x1=1100.0;
    x2=2300.0;
    dx=dxin;
    
    z1=0;
    z2=200;
    dz=dzin;
    
    xs=x1:dx:x2;
    zs=z1:dz:z2;
    
    nDep=length(zs);
    nOff=length(xs);
    
    InputDataParams.skipSta=skipSta;
    InputDataParams.KTimes=KTimes;
    InputDataParams.incangs=incangs;
    InputDataParams.TakeDifferences=TakeDifferences;
    InputDataParams.DeconvolveParentWaveform=DeconvolveParentWaveform;
    
    
    fprintf('Building input matrices... ')
    if direction==1;
        [Data,Locations,RayParams,BackAzimuths] = build_input_matrices(InputDataParams);
    elseif direction == 2;
        [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_reverse(InputDataParams);
    else
        [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_both_directions(InputDataParams);
    end
    fprintf('done.\n')
    
    [Data,Locations,RayParams,BackAzimuths] = random_sample(Data,Locations,RayParams,BackAzimuths,0.0);
    
    nSeis=length(RayParams);
    
    G=zeros(nTimes*nSeis,nOff,nDep);
    %d=zeros(nTimes,1);
    
    %Here we make a matrix C nOff by nDep which store icol values
    C=zeros(nOff,nDep);
    icol=0;
    for ioff = 1:1:nOff;
       for idep = 1:1:nDep;
            icol=icol+1;
            C(ioff,idep)=icol;
       end
    end

    %R is the regularization matrix
    R=sparse(nOff*nDep,nOff*nDep);
    %[X,Y,Z] = meshgrid(-stalocs+1500,Scat_Depths,1:nTimes);
    %volume_input=zeros(nDep,nOff);
    fprintf('Generating G matrix:\n');
    for ii = 1:nSeis;
        Location=Locations(ii);
        BackAzimuth=BackAzimuths(ii);
        RayParam=RayParams(ii);
        %kluge for testing
        if RayParam == 20;
           iRP=1;
        elseif RayParam == 23;
            iRP=2;
        elseif RayParam == 26;
            iRP=3;
        else
            fprintf('***Warning: Bad RP %f \n', RayParam)
        end
        
        fprintf('  %.2f pct complete... \n',(ii-1)/nSeis*100.0);
                
        [Xq,Yq,Zq]=meshgrid(xs-Location,zs,1:nTimes);
        if BackAzimuth==1;
            %[X,Y,Z] = meshgrid(-stalocs+1500,Scat_Depths,1:nTimes);
            tmp=interp3(X,Y,Z,Kernel(:,:,1:nTimes,iRP),Xq,Yq,Zq,'linear',0.0);
        else
            %[X,Y,Z] = meshgrid(+stalocs-1500,Scat_Depths,1:nTimes);
            tmp=interp3(-X,Y,Z,Kernel(:,:,1:nTimes,iRP),Xq,Yq,Zq,'linear',0.0);
        end
        %tmp=interp3(X,Z,Y,Kernel(:,1:nTimes,:),Xq,Zq,Yq,'spline',0.0);
        G((1:nTimes)+(ii-1)*nTimes,:,:)=permute(tmp,[3 2 1]); 
    
        %loop over offsets and depths
        for ioff = 1:1:nOff;
            for idep = 1:1:nDep;

                icol=C(ioff,idep);
                %ioff_for_sta=ioff-(Stas(ista)+1);
                
                %Construct regularization matrix
                if idep > 1 && idep < nDep && ioff > 1 && ioff < nOff;        
                    %icolu=C(ioff,idep+1);
                    %icold=C(ioff,idep-1);
                    icoll=C(ioff+1,idep);
                    icolr=C(ioff-1,idep);
                    R(icol,icol)  = 1.0;
                    %R(icol,icolu) =-0.25;
                    %R(icol,icold) =-0.25;
                    R(icol,icoll) =-0.5;
                    R(icol,icolr) =-0.5;
                elseif idep ==1 && ioff > 1 && ioff < nOff;
                    %icolu=C(ioff,idep+1);
                    icoll=C(ioff+1,idep);
                    icolr=C(ioff-1,idep);
                    R(icol,icol)  = 1.0;
                    %R(icol,icolu) =-1.0/3.0;
                    R(icol,icoll) =-1.0/2.0;
                    R(icol,icolr) =-1.0/2.0;
                elseif idep ==nDep && ioff > 1 && ioff < nOff;
                    %icold=C(ioff,idep-1);
                    icoll=C(ioff+1,idep);
                    icolr=C(ioff-1,idep);
                    R(icol,icol)  = 1.0;
                    %R(icol,icold) =-1.0/3.0;
                    R(icol,icoll) =-1.0/2.0;
                    R(icol,icolr) =-1.0/2.0;
                elseif idep > 1 && idep < nDep && ioff == 1;
                    %icolu=C(ioff,idep+1);
                    %icold=C(ioff,idep-1);
                    icoll=C(ioff+1,idep);
                    R(icol,icol)  = 1.0;
                    %R(icol,icolu) =-1.0/3.0;
                    %R(icol,icold) =-1.0/3.0;
                    R(icol,icoll) =-1.0/1.0;
                elseif idep > 1 && idep < nDep && ioff == nOff;
                    %icolu=C(ioff,idep+1);
                    %icold=C(ioff,idep-1);
                    icolr=C(ioff-1,idep);
                    R(icol,icol)  = 1.0;
                    %R(icol,icolu) =-1.0/3.0;
                    %R(icol,icold) =-1.0/3.0;
                    R(icol,icolr) =-1.0/1.0;
                else
                    %Don't do anything to the corners for now...
                    %0.0 for nothing, 1.0 for small
                    R(icol,icol) = 1.0;
                end
           end            
        end
    end
    
    G=permute(G,[1,3,2]);
    G=reshape(G,nTimes*nSeis,nOff*nDep);
    G=sparse(G);

    ireg=0;
    
    d=reshape(Data',[],1);
    
    
    for norm_opt = Norm_Opts;
        
        if norm_opt==1;
           disp('Minimizing model norm');Rk=sparse(diag(ones(nOff*nDep,1)));
        elseif norm_opt==2;
           Rk=R;
           disp('Minimizing model roughness');
        else
            Rk=R*0;
            disp('Not regularized')
        end
    
        for nu = Nus;
            ireg=ireg+1;

            [dhat,volume,vred]=invert_iterative(G,Rk,nu,d,nOff,nDep,nIterMax);
          
            tmp=sprintf('_%d.mat',ireg);
            filename=[saveFilename tmp];
            
            if (TakeDifferences)
                model.vs=[3.2,4.4,4.05];
                model.hs=[60.0,60.0,180.0];
                for ii = 1:length(xs)
                    for jj = 1:length(zs)
                        model.vs=[3.2,4.4,4.05];
                        dv=get_velocity_from_1Dmod(model,zs(jj));
                        disp(dv);
                        volume(jj,ii)=volume(jj,ii);
                    end
                end
            end
            
            save(filename,'dhat','volume','vred','xs','zs','nu','norm_opt','Locations','nSeis','KTimes','d','Kernel_Type')
            
        end
    end
end

function [Data,Locations,RayParams,BackAzimuths] = build_input_matrices(InputDataParams)

    skipSta=InputDataParams.skipSta;
    KTimes= InputDataParams.KTimes;
    incangs=InputDataParams.incangs;
    TakeDifferences=InputDataParams.TakeDifferences;
    DeconvolveParentWaveform=InputDataParams.DeconvolveParentWaveform;
    

    Sta1=1;
    Sta2=150;
    Stas=Sta1:skipSta:Sta2;
    
    if (skipSta < 0);
        Stas=Sta2:skipSta:Sta1;
    end
    
    nSeis=length(Stas)*length(incangs);
    
    %initialize matrices
    Data=zeros(nSeis,length(KTimes));
    Locations=zeros(nSeis,1);
    RayParams=zeros(nSeis,1);
    BackAzimuths=zeros(nSeis,1);
    
    iseis=0;
    
    
    for incang = incangs;
    
        ModelDirectory=sprintf('TEST_MODELS/L2/OUTPUT_FILES_%d-5000-400000/',incang);
        RfnceDirectory=sprintf('TEST_MODELS/L2/OUTPUT_FILES_%d-0-400000/',incang);
        
        command=['cat ' ModelDirectory '/DATA/STATIONS | awk {''print $3''} > tmp.txt'];
        system(command);
        AllLocations=load('tmp.txt');
        AllLocations=AllLocations/double(1000.0); % convert to km


        for iSta = Stas;
            iseis=iseis+1;
            Daughter=load(sprintf([ModelDirectory 'OUTPUT_FILES/AA.S0%03d.BXP.semd'],iSta));
            Parent  =load(sprintf([ModelDirectory 'OUTPUT_FILES/AA.S0%03d.BXS.semd'],iSta));
            
            DRef    =load(sprintf([RfnceDirectory 'OUTPUT_FILES/AA.S0%03d.BXP.semd'],iSta));
            %PRef    =load(sprintf([RfnceDirectory 'OUTPUT_FILES/AA.S0%03d.BXS.semd'],iSta));
            
            [AmpMax,imax]=max(abs(Parent(:,2)));
            time=Parent(:,1) - Parent(imax,1);

            %Parent(:,2)=Parent(:,2)/AmpMax;
            Daughter(:,2)=Daughter(:,2)/AmpMax;
            
            %plot(Parent(:,1),Parent(:,2),Daughter(:,1),Daughter(:,2))
            %pause
            
            %Subtract reference waveform
            if (TakeDifferences)
                Daughter(:,2)=Daughter(:,2)-DRef(:,2);
                %Parent(:,2)  =Parent(:,2)  -PRef(:,2);
            end
            
            if (DeconvolveParentWaveform)
               Mask=(Parent(:,1)>Parent(imax-40,1)).*(Parent(:,1)<Parent(imax+40,1));
               P=Parent(:,2).*Mask;
               D=Daughter(:,2);
               TB=4;
               NT=3;
               dt=Parent(2,1)-Parent(1,1);
               win_len=50;
               Poverlap=0.99;        
                
               [Time, RF_Time] = ETMTM(P(1:imax+200)',D(1:imax+200)',TB,NT,'data',dt,win_len,Poverlap);
               
               %clf;
               %subplot(2,1,1)
               %plot(P); hold on
               %plot(D);
               %subplot(2,1,2)
               %plot(Time,RF_Time);
               %pause
           
               %reset and resize daughter
               Daughter=zeros(length(RF_Time),2);
               Daughter(:,2)=RF_Time;
               time=Time;
               
            end
            
            dtmp=interp1(time,Daughter(:,2),KTimes)';

            Data(iseis,:)=dtmp;
            BackAzimuths(iseis,1)=1.0;
            RayParams(iseis,1)=incang;
            Locations(iseis,1)=AllLocations(iSta);

        end
    end
    
end

function [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_reverse(InputDataParams)

    InputDataParams.skipSta=-InputDataParams.skipSta;
    
    [Data,Locations,RayParams,BackAzimuths]=build_input_matrices(InputDataParams);
    
    %Reflect locations about symmetry axis
    x0=1725.0;
    tmp= 2*x0 - Locations;
    
    Locations = tmp;
    
    BackAzimuths=BackAzimuths*-1.0;
    
end

function [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_both_directions(InputDataParams)

[Data1,Locations1,RayParams1,BackAzimuths1] = build_input_matrices(InputDataParams);
[Data2,Locations2,RayParams2,BackAzimuths2] = build_input_matrices_reverse(InputDataParams);

Data = [Data1; Data2];
Locations = [Locations1; Locations2];
RayParams = [RayParams1; RayParams2];
BackAzimuths = [BackAzimuths1; BackAzimuths2];

end

function [x,R] = cgstep(x,g,R,G,first)
persistent s
persistent S
if first; %do steepest descents
    s=x*0.0;
    S=R*0.0;
    beta=  0.0;
    alfa= - dot(G,R) / dot(G,G);
else
    gdg=dot(G,G);
    sds=dot(S,S);
    gds=dot(G,S);
    
    determ  = gdg * sds * max(1.0-(gds/gdg)*(gds/sds),1.0E-12);
    
    gdr = - dot(G,R);
    sdr = - dot(S,R);
    alfa = ( sds*gdr - gds*sdr) / determ;
    beta = (-gds*gdr + gdg*sdr) / determ;
end
s  = alfa*g + beta*s;
S  = alfa*G + beta*S;
x  = x + s;
R  = R + S;

end

function [dhat,volume,vred]=invert_iterative(G,R,nu,d,nOff,nDep,nIterMax)
    
    %Full system
    Gfull=(horzcat(G',R*nu)');
    dfull=[d',zeros(1,nOff*nDep)]';
    
    %Back-project
    GT=Gfull';
    m=GT*dfull;
    m=m*0.0;
 
    %Iterate (steepest descents)
    %
    %
    r=Gfull*m-dfull;
    dm=GT*r;
    dr=Gfull*dm;
    g=GT*r;
    [m,r]= cgstep(m,g,r,dr,true);
    
    for iter=1:nIterMax;
        dm=GT*r;
        dr=Gfull*dm;
        g=GT*r;
        [m,r]= cgstep(m,g,r,dr,false);
        vred=1-norm(r(1:length(d)))/norm(d);
        fprintf('Iteration %4d, Var Red = %.4f pct\n',iter,vred*100.0)
        if exist('vred_last','var')
            if vred-vred_last<0.000001;
                break
            end
        end
        vred_last=vred;
    end
    
    dhat=G*m;
    
    icol=0;
    volume=zeros(nDep,nOff);
    for ioff = 1:1:nOff;
        for idep = 1:1:nDep;
           icol=icol+1;
            volume(idep,ioff)=m(icol); 
        end
    end
    
end

function [Data,Locations,RayParams,BackAzimuths] = random_sample(Data,Locations,RayParams,BackAzimuths,factor)

nData=length(Locations);

mask=rand(nData,1)>factor;

    Data=Data(mask,:);
    Locations=Locations(mask);
    RayParams=RayParams(mask);
    BackAzimuths=BackAzimuths(mask);

end
