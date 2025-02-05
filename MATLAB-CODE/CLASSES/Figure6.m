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
                      'OUTPUT/24-Aug-2017/LOOP-%d-%d-%d-%d-%s',...
                      labprops.amplitude,labprops.wavlen,labprops.depth,...
                      1,'false');
                  
                try
                    Inversion=ReadFromDisk(obj);
                    %fprintf('min, max G = %f, %f\n', min(min(Inversion.DataVector.d)), max(max(Inversion.DataVector.d)));
                    %print(Inversion.DataVector)
                    %figure(99);
                    %plot(extractSeis(Inversion.DataVector,75)); hold on;
                    figure(1);
                    obj.fig=plot_model(Inversion.VelocityModel2D(),'figure6',obj.clabel,labprops,ax);
                    caxis([-500, 500])
                    
                    if iplt == 6
                        c=colorbar;

                        totwid=ax.Position(3);

                        gr=1.618;
                        cwid=totwid*0.03;
                        buf=cwid*gr;

                        ax.Position(3)=totwid*1.00;
                        c.Position(1)=ax.Position(1)+totwid+buf;
                        c.Position(2)=ax.Position(2);
                        c.Position(3)=cwid;
                        c.Position(4)=ax.Position(4)/gr;
                        
                        c.Label.String='Scattering Potential';
                    end
                    
                    
                    if mod(iplt,2)==0;
                        ylabel('');
                    end
                    if iplt < 5
                        xlabel('');
                    end
                    title('')
                    
                catch
                    title('No file found')
                end
            end
            suptitle('~5 km Station Spacing')
            save(obj,'Figure6')
        end

    end
    methods (Access = private)        
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
