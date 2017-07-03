classdef isochrons
    properties
        Ps
        Sp
        xs
        zs
        rp
        model
    end
    methods
        function obj=isochrons(xs,zs,rp,model)
            obj.xs=xs;
            obj.zs=zs;
            obj.rp=rp;
            obj.model.vp=model.vp;
            obj.model.vs=4.4;
            obj.model.hs=model.hs;
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