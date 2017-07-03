classdef isochrons
    properties
        Sp
        xs
        zs
        rp
    end
    methods
        function obj=isochrons(xs,zs,rp)
            obj.xs=xs;
            obj.zs=zs;
            obj.rp=rp;
            model.vp=5.0;
            model.vs=3.0;
            model.hs=300.0;
         
            [Tdirect,Tscat,~,~,~]=read_in_necessary_files(rp, xs, zs, model);
            [T] = timeshifts(Tdirect,Tscat,xs,zs);
         
            obj.Sp=T;
        end 
    end
end