classdef Figure6 < MyFigure
    properties
        clabel='Scattering Potential'
        path='OUTPUT/07-Jul-2017/1'
    end
    methods
        function obj=Figure6()
            Inversion=ReadFromDisk(obj);
            plot_model(Inversion.VelocityModel2D(),'figure6',obj.clabel)
        end

    end
    methods (Access = private)
        
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
            
    end
    
end
