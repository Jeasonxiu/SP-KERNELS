classdef Figure7 < MyFigure 
    properties
        clabel='Scattering Potential'
        path
    end
    methods
        function obj=Figure7()
            figure(1); clf;
            labprops.wavlen=400000;
            labprops.depth=120000;
            labprops.amplitude=5000;
            
            for iplt = 1:6;
                ax=subplot(3,2,iplt);
                if iplt == 1
                    skipsta=1;
                    GaussianFilter='false';
                elseif iplt == 3
                    skipsta=2;
                    GaussianFilter='false';
                elseif iplt == 5
                    skipsta=4;
                    GaussianFilter='false';
                elseif iplt == 2
                    skipsta=1;
                    GaussianFilter='true';
                elseif iplt == 4
                    skipsta=2;
                    GaussianFilter='true';
                elseif iplt == 6
                    skipsta=4;
                    GaussianFilter='true';
                end
                    obj.path=sprintf(...
                      'OUTPUT/24-Aug-2017/LOOP-%d-%d-%d-%d-%s',...
                      labprops.amplitude,labprops.wavlen,labprops.depth,...
                      skipsta,GaussianFilter);
                  
                try
                    Inversion=ReadFromDisk(obj);
                    obj.fig=plot_model(Inversion.VelocityModel2D(),'figure6',obj.clabel,labprops,ax);
                    
                    if iplt >= 5
                        c=colorbar;

                        totwid=ax.Position(3);

                        gr=1.618;
                        cwid=totwid*0.03;
                        buf=cwid*gr;

                        c.Position(1)=ax.Position(1)+totwid+buf;
                        c.Position(2)=ax.Position(2);
                        c.Position(3)=cwid;
                        c.Position(4)=ax.Position(4)/gr;
                        

                    end
                    
                    if iplt == 6
                        c.Label.String='Scattering Potential';
                    end
                        
                    
                    
                    if mod(iplt,2)==0;
                        ylabel('');
                    end
                    if iplt < 5
                        xlabel('');
                    end

                    
                catch
                    fprintf('%s not found\n',obj.path);
                    title('No file found')
                end
                
                
            end
            save(obj,'Figure7')
        end

    end
    methods (Access = private)        
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
