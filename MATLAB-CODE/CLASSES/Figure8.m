classdef Figure8 < MyFigure 
    properties
        clabel='Scattering Potential'
        path
    end
    methods
        function obj=Figure8()
            fig=figure(1); clf;
            fig.Units='Inches';
            fig.Position(1)=1;
            fig.Position(2)=2;
            fig.Position(3)=4.5;
            fig.Position(4)=8;
            for iplt = 1:3;
                ax=subplot(3,1,iplt);
                if iplt == 1
                    labprops.wavlen=200000;
                    labprops.depth=120000;
                    labprops.amplitude=10000;
                elseif iplt == 2
                    labprops.wavlen=400000;
                    labprops.depth=120000;
                    labprops.amplitude=20000;
                elseif iplt == 3
                    labprops.wavlen=100000;
                    labprops.depth=120000;
                    labprops.amplitude=10000;
                end
                    obj.path=sprintf(...
                      'OUTPUT/24-Aug-2017/LOOP-%d-%d-%d-%d-%s',...
                      labprops.amplitude,labprops.wavlen,labprops.depth,...
                      1,'false');
                  
                try
                    Inversion=ReadFromDisk(obj);
                    obj.fig=plot_model(Inversion.VelocityModel2D(),'figure8',obj.clabel,labprops,ax);
                    caxis([-500, 500])
                    
                    ax.Position(3)=ax.Position(3)*0.90;
                    
                    
                    if iplt == 3
                        totwid=ax.Position(3);
                        
                        c=colorbar;



                        gr=1.618;
                        cwid=totwid*0.03;
                        buf=cwid*gr;

                        c.Position(1)=ax.Position(1)+totwid+buf;
                        c.Position(2)=ax.Position(2)+buf;
                        c.Position(3)=cwid;
                        c.Position(4)=ax.Position(4)/gr;
                        c.Label.String='Scattering Potential';
                    end
                    
                    if iplt == 2 || iplt == 3;
                        title('');
                    end
                    
                    if iplt == 1 || iplt == 2;
                        xlabel('');
                    end
                    
                catch
                    title('No file found')
                end
            end
            save(obj,'Figure8')
        end

    end
    methods (Access = private)        
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
