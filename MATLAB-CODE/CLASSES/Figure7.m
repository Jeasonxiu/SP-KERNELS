classdef Figure7 < MyFigure
    properties
        clabel='Velocity Perturbation (%)'
        path='OUTPUT/28-Jul-2017/NEWNORM2'
    end
    methods
        function obj=Figure7()
            Inversion=ReadFromDisk(obj);
            plot_model(Inversion.VelocityModel2D(),'figure7',obj.clabel)
        end

    end
    methods (Access = private)
        
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
