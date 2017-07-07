classdef Inversion
    properties
        InversionParams=InversionParams();
        velocity_model=velocity_model();
        kernel=kernel(3);
        VelocityModel2D=VelocityModel2D();
        d=[]
        nOff=[]
        nDep=[]
    end
    properties (Hidden)
        G
        R
    end
    methods 
        function obj=Inversion(varargin)
            if (nargin>0)
                obj=SetDefaultInversionParams(obj);
                obj=SetUpKernel(obj);
                obj=SetUpMatrices(obj);
                obj=RunInversion(obj);
                plot_model(obj.VelocityModel2D)
    
            end
        end
        function obj=SetDefaultInversionParams(obj)
            obj.InversionParams=SetDefaultParams2(obj.InversionParams);
        end
        function obj=SetUpKernel(obj)
            if (obj.InversionParams.ImagingMethod == 1)
                fprintf('Performing back-projection migration\n')
                %Thin kernels
                obj.kernel.tchar=0.1;
                obj.kernel.nderiv=0.0;
                %One iteration
                obj.InversionParams.nIterMax=1;
                %Don't take data difference
                obj.InversionParams.TakeDifferences=false;
                %Don't regularize
                obj.InversionParams.Norm_Opts=0;
                %Don't Deconvolve
                obj.InversionParams.DeconvolveParentWaveform=false;

            elseif (obj.InversionParams.ImagingMethod == 2)
                fprintf('Performing C-G Inversion\n')
                %FF kernels
                obj.kernel.tchar=0.8;
                obj.kernel.nderiv=3.0;
                %Many iterations
                %use obj.InversionParams.nIterMax
                %Take data difference
                obj.InversionParams.TakeDifferences=true;
                %Don't regularize
                obj.InversionParams.Norm_Opts=2;
                %Don't deconvolve
                obj.InversionParams.DeconvolveParentWaveform=false;
            else
                return
            end

            obj.kernel.Angles=[20,23,26];
            
            disp(obj.kernel.Angles)
            
            obj.kernel=load(obj.kernel);
            
            %Kernel.kernel=real(Kernel);
            %Kernel.kernel(isnan(Kernel))=0;
            
        end  
        function obj=SetUpMatrices(obj)
            
            %Subfunctions
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
            
            function [Data,Locations,RayParams,BackAzimuths] = random_sample(Data,Locations,RayParams,BackAzimuths,factor)

                nData=length(Locations);

                mask=rand(nData,1)>factor;

                Data=Data(mask,:);
                Locations=Locations(mask);
                RayParams=RayParams(mask);
                BackAzimuths=BackAzimuths(mask);

            end
            
            %Set model parameters
            x1=1100.0;
            x2=2300.0;
            dx=obj.InversionParams.dxin;
    
            z1=0;
            z2=200;
            dz=obj.InversionParams.dzin;
    
            obj.VelocityModel2D.xs=x1:dx:x2;
            obj.VelocityModel2D.zs=z1:dz:z2;

            obj.nDep=length(obj.VelocityModel2D.zs);
            obj.nOff=length(obj.VelocityModel2D.xs);

            InputDataParams.skipSta=obj.InversionParams.skipSta;
            InputDataParams.KTimes=obj.kernel.KTimes;
            InputDataParams.incangs=obj.kernel.Angles;
            InputDataParams.TakeDifferences=obj.InversionParams.TakeDifferences;
            InputDataParams.DeconvolveParentWaveform=obj.InversionParams.DeconvolveParentWaveform;


            fprintf('Building input matrices... ')
            if obj.InversionParams.direction==1;
                [Data,Locations,RayParams,BackAzimuths] = build_input_matrices(InputDataParams);
            elseif obj.InversionParams.direction == 2;
                [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_reverse(InputDataParams);
            else
                [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_both_directions(InputDataParams);
            end
            fprintf('done.\n')

            [Data,Locations,RayParams,BackAzimuths] = random_sample(Data,Locations,RayParams,BackAzimuths,0.0);

            obj.VelocityModel2D.Locations=Locations;
            
            nSeis=length(RayParams);
            
            obj.VelocityModel2D.nSeis=nSeis;
            
            obj.G=zeros(obj.kernel.nTimes*nSeis,obj.nOff,obj.nDep);
            %d=zeros(nTimes,1);

            %Here we make a matrix C obj.nOff by obj.nDep which store icol values
            C=zeros(obj.nOff,obj.nDep);
            icol=0;
            for ioff = 1:1:obj.nOff;
               for idep = 1:1:obj.nDep;
                    icol=icol+1;
                    C(ioff,idep)=icol;
               end
            end

            %R is the regularization matrix
            obj.R=sparse(obj.nOff*obj.nDep,obj.nOff*obj.nDep);
            %[X,Y,Z] = meshgrid(-stalocs+1500,Scat_Depths,1:nTimes);
            %volume_input=zeros(obj.nDep,obj.nOff);
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

                [Xq,Yq,Zq]=meshgrid(obj.VelocityModel2D.xs-Location,obj.VelocityModel2D.zs,1:obj.kernel.nTimes);
                [obj.kernel.X,obj.kernel.Y,obj.kernel.Z]=meshgrid(obj.kernel.xs,obj.kernel.zs,1:obj.kernel.nTimes);
                
                if BackAzimuth==1;
                    %[X,Y,Z] = meshgrid(-stalocs+1500,Scat_Depths,1:nTimes);
                    tmp=interp3(obj.kernel.X,obj.kernel.Y,obj.kernel.Z,obj.kernel.Kernel(:,:,1:obj.kernel.nTimes,iRP),Xq,Yq,Zq,'linear',0.0);
                else
                    %[X,Y,Z] = meshgrid(+stalocs-1500,Scat_Depths,1:nTimes);
                    tmp=interp3(-obj.kernel.X,obj.kernel.Y,obj.kernel.Z,obj.kernel.Kernel(:,:,1:obj.kernel.nTimes,iRP),Xq,Yq,Zq,'linear',0.0);
                end
                %tmp=interp3(X,Z,Y,Kernel(:,1:nTimes,:),Xq,Zq,Yq,'spline',0.0);
                nTimes=obj.kernel.nTimes;
                obj.G((1:nTimes)+(ii-1)*nTimes,:,:)=permute(tmp,[3 2 1]); 
                
                
                
                %loop over offsets and depths
                for ioff = 1:1:obj.nOff;
                    for idep = 1:1:obj.nDep;

                        icol=C(ioff,idep);
                        %ioff_for_sta=ioff-(Stas(ista)+1);

                        %Construct regularization matrix
                        if idep > 1 && idep < obj.nDep && ioff > 1 && ioff < obj.nOff;        
                            %icolu=C(ioff,idep+1);
                            %icold=C(ioff,idep-1);
                            icoll=C(ioff+1,idep);
                            icolr=C(ioff-1,idep);
                            obj.R(icol,icol)  = 1.0;
                            %obj.R(icol,icolu) =-0.25;
                            %obj.R(icol,icold) =-0.25;
                            obj.R(icol,icoll) =-0.5;
                            obj.R(icol,icolr) =-0.5;
                        elseif idep ==1 && ioff > 1 && ioff < obj.nOff;
                            %icolu=C(ioff,idep+1);
                            icoll=C(ioff+1,idep);
                            icolr=C(ioff-1,idep);
                            obj.R(icol,icol)  = 1.0;
                            %obj.R(icol,icolu) =-1.0/3.0;
                            obj.R(icol,icoll) =-1.0/2.0;
                            obj.R(icol,icolr) =-1.0/2.0;
                        elseif idep ==obj.nDep && ioff > 1 && ioff < obj.nOff;
                            %icold=C(ioff,idep-1);
                            icoll=C(ioff+1,idep);
                            icolr=C(ioff-1,idep);
                            obj.R(icol,icol)  = 1.0;
                            %obj.R(icol,icold) =-1.0/3.0;
                            obj.R(icol,icoll) =-1.0/2.0;
                            obj.R(icol,icolr) =-1.0/2.0;
                        elseif idep > 1 && idep < obj.nDep && ioff == 1;
                            %icolu=C(ioff,idep+1);
                            %icold=C(ioff,idep-1);
                            icoll=C(ioff+1,idep);
                            obj.R(icol,icol)  = 1.0;
                            %obj.R(icol,icolu) =-1.0/3.0;
                            %obj.R(icol,icold) =-1.0/3.0;
                            obj.R(icol,icoll) =-1.0/1.0;
                        elseif idep > 1 && idep < obj.nDep && ioff == obj.nOff;
                            %icolu=C(ioff,idep+1);
                            %icold=C(ioff,idep-1);
                            icolr=C(ioff-1,idep);
                            obj.R(icol,icol)  = 1.0;
                            %obj.R(icol,icolu) =-1.0/3.0;
                            %obj.R(icol,icold) =-1.0/3.0;
                            obj.R(icol,icolr) =-1.0/1.0;
                        else
                            %Don't do anything to the corners for now...
                            %0.0 for nothing, 1.0 for small
                            obj.R(icol,icol) = 1.0;
                        end
                   end            
                end
            end

            obj.G=permute(obj.G,[1,3,2]);
            obj.G=reshape(obj.G,obj.kernel.nTimes*nSeis,obj.nOff*obj.nDep);
            obj.G=sparse(obj.G);

            

            obj.d=reshape(Data',[],1);
        end
        function obj=RunInversion(obj)
            %subfunctions
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
            
            ireg=0;
            
            for norm_opt = obj.InversionParams.Norm_Opts;
        
                if norm_opt==1;
                   disp('Minimizing model norm');Rk=sparse(diag(ones(obj.nOff*obj.nDep,1)));
                elseif norm_opt==2;
                   Rk=obj.R;
                   disp('Minimizing model roughness');
                else
                    Rk=obj.R*0;
                    disp('Not regularized')
                end

                for nu = obj.InversionParams.Nus;
                    ireg=ireg+1;

                    [dhat,volume,vred]=invert_iterative(obj.G,Rk,nu,obj.d,obj.nOff,obj.nDep,obj.InversionParams.nIterMax);

                    if (obj.InversionParams.TakeDifferences)
                        for ii = 1:length(obj.VelocityModel2D.xs)
                            for jj = 1:length(obj.VelocityModel2D.zs)
                                [~,dv]=get_v(obj.velocity_model,obj.VelocityModel2D.zs(jj));
                                volume(jj,ii)=volume(jj,ii);
                            end
                        end
                    end

                    %save(filename,'dhat','volume','vred','xs','zs','nu','norm_opt','Locations','nSeis','KTimes','d','Kernel_Type')

                    obj.VelocityModel2D.dlnvs=volume;
                    obj.VelocityModel2D.dhat=dhat;
                    obj.VelocityModel2D.vred=vred;
                    obj.VelocityModel2D.nu=nu;
                    obj.VelocityModel2D.d=obj.d;
                    obj.VelocityModel2D.dtime=obj.kernel.KTimes;
                    
                end
            end
            
            
        end
        function obj=PurgeLargeMatrices(obj)
            obj.G=[];
            obj.R=[];
        end
        function Save(obj,varargin)
            if nargin>1
                filename=varargin{1};
            else
                filename='I';
            end
            %G and R too large to save
            Inversion=PurgeLargeMatrices(obj); %#ok<NASGU>
            save(filename,'Inversion')
        end
    end
    
end