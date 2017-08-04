classdef Inversion
    properties
        InversionParams=InversionParams();
        velocity_model=velocity_model();
        kernel=kernel(3);
        VelocityModel2D=VelocityModel2D();
        DataParams=[]
        DataVector=[]
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
                obj=SetDefaultInversionParams(obj, 1);
                obj=SetUpKernel(obj);
                obj=SetUpMatrices(obj);
                obj=RunInversion(obj);
                plot_model(obj.VelocityModel2D)
            end
        end
        function obj=SetDefaultInversionParams(obj,iParam)
            if iParam == 1
                obj.InversionParams=SetDefaultParams1(obj.InversionParams);
            elseif iParam == 2
                obj.InversionParams=SetDefaultParams2(obj.InversionParams);
            else
                error('Bad iParam: %d', iParams)
            end
        end
        function obj=SetUpKernel(obj)
            if (obj.InversionParams.ImagingMethod == 1)
                fprintf('Performing back-projection migration\n')
                %Thin kernels
                obj.kernel.tchar=0.20;
                obj.kernel.nderiv=0.5;
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
        function obj=SetDataParams(obj,lab_amplitude,lab_wavelength,lab_depth,skipSta)
            obj.DataParams=DataParams(lab_amplitude,lab_wavelength,lab_depth,...
                skipSta,obj.kernel.KTimes,...
                obj.kernel.Angles,obj.InversionParams.TakeDifferences,...
                obj.InversionParams.DeconvolveParentWaveform,...
                obj.InversionParams.direction);
        end
        function obj=SetUpMatrices(obj)
            
            DV=DataVector(obj.DataParams);
            obj.DataVector=DV;
            
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

            obj.VelocityModel2D.Locations=DV.Locations;
            
            nSeis=length(DV.RayParams);
            
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
                Location=DV.Locations(ii);
                BackAzimuth=DV.BackAzimuths(ii);
                RayParam=DV.RayParams(ii);
                %kluge for testing
                for iRP = 1:length(DV.DataParams.Angles);
                    if RayParam == DV.DataParams.Angles(iRP);
                        break
                    end
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
                
                volume=obj.m2volume(m,obj);

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

                    [dhat,volume,vred]=invert_iterative(obj.G,Rk,nu,obj.DataVector.d,obj.nOff,obj.nDep,obj.InversionParams.nIterMax);

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
                    obj.VelocityModel2D.d=obj.DataVector.d;
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
    
    methods (Static)       
        function volume=m2volume(m,obj)
            icol=0;
            volume=zeros(obj.nDep,obj.nOff);
            for ioff = 1:1:obj.nOff;
                for idep = 1:1:obj.nDep;
                   icol=icol+1;
                   volume(idep,ioff)=m(icol); 
                end
            end
        end
        function m=volume2m(volume,obj)
                icol=0;
                m=zeros(obj.nDep*obj.nOff,1);
                for ioff = 1:1:obj.nOff;
                    for idep = 1:1:obj.nDep;
                       icol=icol+1;
                       m(icol)=volume(idep,ioff); 
                    end
                end

        end
        
        function obj=testForwardModel(obj)
            volume=obj.generateTestVolume(obj);
            m=obj.volume2m(volume,obj);
            dhat=obj.G*m;
         
            
            obj.VelocityModel2D.dlnvs=volume;
            obj.VelocityModel2D.dhat=dhat;
            
        end
        
        function volume=generateTestVolume(obj)
            dsmooth=1; %km
            volume=zeros(obj.nDep,obj.nOff);
            for ioff = 1:1:obj.nOff;
                for idep = 1:1:obj.nDep;
                   distFromLayer=60-obj.VelocityModel2D.zs(idep);
                   
                   fac1=exp(-distFromLayer^2/2/dsmooth^2);
                   fac2=(1-(distFromLayer/dsmooth)^2);
                   amp=-3000*fac1;
                   
                   volume(idep,ioff)=amp; 
                end
            end 
        end
        
    end
    
end
