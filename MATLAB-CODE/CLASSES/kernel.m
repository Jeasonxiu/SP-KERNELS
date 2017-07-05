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
            obj.xs=(-600:5:150)';
            obj.zs=0:5:300;
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
        function [G,XSP,A]=get_amplitude_components(obj,Pdirect,timeref)
            if obj.Kernel_Type~=3;
                error('Bad Kernel_Type: %d',obj.Kernel_Type);
            else
                [Tdirect,Tscat,Pscat,~,~]=read_in_necessary_files(Pdirect, obj.xs, obj.zs, obj.model);
                [G] = geom_spreading(Pscat,obj.xs,obj.zs,obj.model);
                [THE,~,~] = calc_angle(Pdirect,Pscat,obj.xs,obj.zs,obj.model);
            %reference speeds 
                vp=7.920;
                vs=4.400;
                rho=3.300;
                vpert=0.01;
                XSP = scattering_pattern(THE,vp,vs,rho,vpert);

                trial_times=-10:0.1:10;
                amps=zeros(length(trial_times),1);
        
                %calculate reference amplitude
                for itime=1:length(trial_times);
                   amps(itime)=fractional_deriv(obj.tchar,2.0,trial_times(itime)); 
                end

                AmpRef=max(abs(amps));

                %calculate actual amplitude
                for itime=1:length(trial_times);
                   amps(itime)=fractional_deriv(obj.tchar,obj.nderiv,trial_times(itime)); 
                end

                amps=amps/AmpRef;

                [T] = timeshifts(Tdirect,Tscat,obj.xs,obj.zs);
                
                time=timeref;  

                %Step 2: Subtract target time from array of times

                DT=T-time;

                %Step 3: Function maps array of times to array of amplitudes (from src time function)
                A = time2amp(DT,trial_times,amps);
                
            end
        end      
    end
end
