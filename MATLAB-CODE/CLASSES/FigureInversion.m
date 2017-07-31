classdef FigureInversion < MyFigure 
    properties
        clabel='Scattering Potential'
        path
        labdepth
        labwavlen
        labamp
    end
    methods
        function obj=FigureInversion(labdepth,labwavlen,labamp)
            obj.labdepth=labdepth;
            obj.labwavlen=labwavlen;
            obj.labamp=labamp;
            obj.path=sprintf(...
                'OUTPUT/31-Jul-2017/LOOP-%d-%d-%d-1',...
                obj.labamp,obj.labwavlen,obj.labdepth);
            Inversion=ReadFromDisk(obj);
            obj.fig=plot_model(...
                Inversion.VelocityModel2D(),'tmp',obj.clabel,...
                obj.labwavlen,obj.labamp,obj.labdepth);
        end

    end
    methods (Access = private)
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
