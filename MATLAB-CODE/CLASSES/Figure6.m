classdef Figure6 < MyFigure 
    properties
        clabel='Scattering Potential'
        path
    end
    methods
        function obj=Figure6()
            figure(1); clf;
            for iplt = 1:6;
                ax=subplot(3,2,iplt);
                if iplt == 1
                    labprops.wavlen=400000;
                    labprops.depth=120000;
                    labprops.amplitude=5000;
                elseif iplt == 3
                    labprops.wavlen=200000;
                    labprops.depth=120000;
                    labprops.amplitude=5000;
                elseif iplt == 5
                    labprops.wavlen=400000;
                    labprops.depth=120000;
                    labprops.amplitude=10000;
                elseif iplt == 2
                    labprops.wavlen=400000;
                    labprops.depth=180000;
                    labprops.amplitude=5000;
                elseif iplt == 4
                    labprops.wavlen=200000;
                    labprops.depth=180000;
                    labprops.amplitude=5000;
                elseif iplt == 6
                    labprops.wavlen=400000;
                    labprops.depth=180000;
                    labprops.amplitude=10000;
                end
                    obj.path=sprintf(...
                      'OUTPUT/17-Aug-2017/LOOP-%d-%d-%d-%d-%s',...
                      labprops.amplitude,labprops.wavlen,labprops.depth,...
                      1,'false');
                  
                try
                    Inversion=ReadFromDisk(obj);
                    obj.fig=plot_model(Inversion.VelocityModel2D(),'figure6',obj.clabel,labprops,ax);
                catch
                    title('No file found')
                end
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
