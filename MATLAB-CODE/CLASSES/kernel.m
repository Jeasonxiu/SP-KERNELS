classdef kernel
    properties
        Kernel
        X
        Y
        Z
        KTimes
        nTimes
        Angles
        tDeci=1;
        Kernel_Type
        tchar=1.0
        nderiv=3.0
    end
    methods
        function obj=kernel(Kernel_Type)
            obj.Kernel_Type=Kernel_Type;
            if Kernel_Type == 1;
                path='/Volumes/nmancine/data2/nmancine/PROJECTS/SP_RECEIVER_FUNCTIONS/KERNEL/SP-KERNELS/DATA';
                [Kernel,X,Y,Z,KTimes,nTimes,Angles] = load_kernel(path,obj.tDeci);
            %elseif Kernel_Type ==2;
            %    [Kernel,zs,xs,~,KTimes] = analytical_kernel_halfspace(26);
            %    nTimes=length(KTimes);
            %    [X,Y,Z]=meshgrid(xs,zs,1:nTimes);
            %    TMP(:,:,:)=Kernel(:,:,1,:);
            %    Kernel=TMP;
            elseif Kernel_Type ==3;
                [~,betatmp]=get_velocity_from_profile('MIGRA/myvmod.nd',300);
                incangs=[20,23,26];
                rps=sind(incangs)/betatmp;
                fprintf('Computing analytical kernel... ')
                model.hs=[60.0,60.0,180.0];
                model.vp=[5.660,7.920,7.280];
                model.vs=[3.2,4.4,4.05];
                xs=(-600:1:150)';
                zs=0:1:300;
                [Kernel,X,Y,Z,KTimes,nTimes] = analytical_kernel_layered(rps, xs, zs, obj.tchar, obj.nderiv, model);
                fprintf('done.\n')
                Angles=incangs;

            else
               error('Bad Kernel_Type: %d',Kernel_Type);
            end
            obj.Kernel=Kernel;
            obj.X=X;
            obj.Y=Y;
            obj.Z=Z;
            obj.KTimes=KTimes;
            obj.nTimes=nTimes;
            obj.Angles=Angles;
        end
    end      
    
end