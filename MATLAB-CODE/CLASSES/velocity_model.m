classdef velocity_model
    %A class containing the velocity model.
    %
    %
    properties
        hs=[60.0,60.0,180.0];
        vp=[5.660,7.920,7.280];
        vs=[3.2,4.4,4.05];
        zmin
        zmax
        nlay;
    end
    methods
        function obj=velocity_model()
            obj=init(obj);
        end
        
        function obj=init(obj)
            obj.nlay=length(obj.hs);
            obj.zmin=0;
            obj.zmax=sum(obj.hs)+obj.zmin;
        end
        
        function obj=update(obj,newmodel)
            obj.hs=newmodel.hs;
            obj.vp=newmodel.vp;
            obj.vs=newmodel.vs;
            obj=init(obj);
        end
        
        function print(obj)
            fprintf('%10s | %10s | %10s \n','H (km)', 'Vp (km/s)', 'Vs (km/s)')
            for ilay=1:obj.nlay;
                fprintf('%10.1f | %10.3f | %10.3f \n',obj.hs(ilay),obj.vp(ilay),obj.vs(ilay))
            end
        end
     
        
        function [vp,vs]=get_v(obj,zs)
            %Get model velocities at depth points
            %
            %
            vp=zeros(1,length(zs))-999.0;
            vs=zeros(1,length(zs))-999.0;     
            
            for iz = 1:length(zs);
                z=zs(iz);
                %check z
                if (z<obj.zmin)
                    error('Too small z= %f', z)
                elseif (z>obj.zmax)
                    error('Too large z= %f', z)    
                end
                
                zref=0;
                for ilay=1:obj.nlay;
                    top=zref;
                    bot=zref+obj.hs(ilay);
                    zref=bot;
                    if ( z>bot )    %depth below layer
                       %
                    elseif (z<top)  %depth above layer
                       %
                    else             %depth within layer
                       vp(iz)=obj.vp(ilay);
                       vs(iz)=obj.vs(ilay);
                       break
                    end

                end

            end

    return

end
    end
end