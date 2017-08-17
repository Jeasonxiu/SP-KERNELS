classdef Figure6 < MyFigure 
    properties
        clabel='Scattering Potential'
        path='OUTPUT/28-Jul-2017/TEST'
    end
    methods
        function obj=Figure6()
            Inversion=ReadFromDisk(obj);
            obj.fig=plot_model(Inversion.VelocityModel2D(),'figure6',obj.clabel);
        end

    end
    methods (Access = private)
        
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
