classdef Figure6 < MyFigure 
    properties
        clabel='Scattering Potential'
        path='OUTPUT/17-Aug-2017/TEST'
    end
    methods
        function obj=Figure6()
            Inversion=ReadFromDisk(obj);
            figure(1); clf;
            for iplt = 1:4;
                ax=subplot(3,2,iplt);
                labprops.wavlen=400000;
                labprops.depth=120000;
                labprops.amplitude=10000;
            	obj.fig=plot_model(Inversion.VelocityModel2D(),'figure6',obj.clabel,labprops,ax);
            end
        end

    end
    methods (Access = private)        
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
