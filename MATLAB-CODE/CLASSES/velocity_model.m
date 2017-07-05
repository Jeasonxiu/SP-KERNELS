classdef velocity_model
    properties
        hs=[60.0,60.0,180.0];
        vp=[5.660,7.920,7.280];
        vs=[3.2,4.4,4.05];
        zmin
        zmax
        nlay;
    end
    methods
        function obj=model()
            obj.nlay=length(obj.hs);
            obj.zmin=0;
            obj.zmax=sum(obj.hs)+obj.zmin;
        end
    end
end