classdef kernel
    properties
        Kernel
        xs
        zs
        rps
        KTimes
        nTimes
        Angles
        Kernel_Type
        model
        X
        Y
        Z
        tDeci=1;
        tchar=1.0
        nderiv=3.0
    end
    methods
        function obj=kernel(Kernel_Type)
            obj.Kernel_Type=Kernel_Type;
            obj.model.hs=[60.0,60.0,180.0];
            obj.model.vp=[5.660,7.920,7.280];
            obj.model.vs=[3.2,4.4,4.05];
            obj.xs=(-600:10:150)';
            obj.zs=0:10:300;
            incangs=[20,23,26];
            [~,betatmp]=get_velocity_from_profile('MIGRA/myvmod.nd',300);
            obj.rps=sind(incangs)/betatmp;
            obj.Angles=incangs;
        end
        
        function obj=load(obj)
            if obj.Kernel_Type == 1;
                path='/Volumes/nmancine/data2/nmancine/PROJECTS/SP_RECEIVER_FUNCTIONS/KERNEL/SP-KERNELS/DATA';
                [obj.Kernel,obj.X,obj.Y,obj.Z,obj.KTimes,obj.nTimes,obj.Angles] = load_kernel(path,obj.tDeci);
            elseif obj.Kernel_Type ==3;
                fprintf('Computing analytical kernel... ')
                [obj.Kernel,~,~,~,obj.KTimes,obj.nTimes] = analytical_kernel_layered(obj.rps, obj.xs, obj.zs, obj.tchar, obj.nderiv, obj.model);
                fprintf('done.\n')
            else
               error('Bad Kernel_Type: %d',obj.Kernel_Type);
            end
            
        end
        
        function [G,XSP]=get_amplitude_components(obj,Pdirect)
            if obj.Kernel_Type~=3;
                error('Bad Kernel_Type: %d',obj.Kernel_Type);
            else
                [~,~,Pscat,~,~]=read_in_necessary_files(Pdirect, obj.xs, obj.zs, obj.model);
                [G] = geom_spreading(Pscat,obj.xs,obj.zs,obj.model);
                [THE,~,~] = calc_angle(Pdirect,Pscat,obj.xs,obj.zs,obj.model);
            %reference speeds 
                vp=7.920;
                vs=4.400;
                rho=3.300;
                vpert=0.01;
                XSP = scattering_pattern(THE,vp,vs,rho,vpert);
            end
        end      
    end
end
