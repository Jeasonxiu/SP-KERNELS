classdef isochrons
    properties
        Ps
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
            model.vp=7.9;
            model.vs=4.4;
            model.hs=300.0;
            [Tdirect,Tscat,~,~,~]=read_in_necessary_files(rp, xs, zs, model);
            [T] = timeshifts(Tdirect,Tscat,xs,zs);
            obj.Sp=T';
            %for Ps, swap velocities
            tmp=model.vp;
            model.vp=model.vs;
            model.vs=tmp;
            [Tdirect,Tscat,~,~,~]=read_in_necessary_files(rp, xs, zs, model);
            [T] = timeshifts(Tdirect,Tscat,xs,zs);
            obj.Ps=T';
            
        end 
    end
end