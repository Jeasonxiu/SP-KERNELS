classdef FigureRecovery < MyFigure 
    properties
        clabel='Scattering Potential'
        path
        labdepth
        labwavlen
        labamp
    end
    methods
        function obj=FigureRecovery()
            
            fig=figure(1); clf;
            fig.Units='Inches';
            fig.Position=[6.0972 7.3889 8.3889 3.6806];
            
            date='24-Aug-2017';
            
            
            syms={'^','o','d'};
            colors={'r','g','b'};
            
            isub=0;
            for labdepth = [120000 180000]
                    isub=isub+1;
                    icurve=0;
                    for labamp = [5000 10000 20000]   
                        icurve=icurve+1;
                        amprecov=nan(3,1);
                        ampstd=nan(3,1);
                        amp10=nan(3,1);
                        amp90=nan(3,1);                   
                        
                        ipoint=0;
                        for labwavlen = [100000 200000 400000]
                            ipoint=ipoint+1;
                        
                            path=sprintf('OUTPUT/%s/LOOP-%d-%d-%d-1-false',date,labamp,labwavlen,labdepth);

                            obj.labdepth=labdepth;
                            obj.labwavlen=labwavlen;
                            obj.labamp=labamp;
                            obj.path=path;
                            labprops.wavlen=obj.labwavlen;
                            labprops.amplitude=obj.labamp;
                            labprops.depth=obj.labdepth;


                            try
                                Inversion=ReadFromDisk(obj);
                                [therat,medrat,stdrat,quant10,quant90] = plot_recovery_performance(Inversion.VelocityModel2D());
                                fprintf('%.3f   %.3f+/-%.3f\n',therat, medrat, stdrat)
                                amprecov(ipoint)=medrat;
                                ampstd(ipoint)=stdrat;
                                amp10(ipoint)=quant10;
                                amp90(ipoint)=quant90;

                            catch
                                fprintf('Error loading.\n')
                            end
                        
                        end
                        
                        subplot(1,2,isub)
                        tmp1=sprintf('-%s%s',syms{icurve},colors{icurve});
                        tmp2=sprintf('-%s',colors{icurve});
                        plot(amprecov,tmp1); hold on
                        %plot(amp90,tmp2); hold on
                        plot([0.8 3.2],[therat therat],'--k','LineWidth',3)
                        subplot(1,2,isub)
                        for itmp = 1:length(ampstd);
                            plot([itmp itmp],[amprecov(itmp)-ampstd(itmp),amprecov(itmp)+ampstd(itmp)],tmp2); hold on
                        end
                        
                        
                    end
                
            end
            
            for isub = 1:2;
                ax=subplot(1,2,isub);
                ax.XTick=[1,2,3];
                ax.XTickLabel=[100,200,400];
                ax.XLim=[0.8,3.2];
            end
            
            for isub = [1 2]
                subplot(1,2,isub);
                xlabel('LAB Wavelength (km)')
            end
            
            for isub = [1 2]
                ax=subplot(1,2,isub);
                ax.YLim=[0, 0.35];
            end
            
            for isub = 1
                subplot(1,2,isub);
                ylabel('Ratio of Amplitude Swings')
            end
            
            %for isub = 3
            %    subplot(2,2,isub)
            %    ylabel('St. Dev. in Ratio of Velocity Swings')
            %end
            
            for isub = 1
                subplot(1,2,isub)
                title('LAB Depth at 120 km')
            end
            
            for isub = 2
                subplot(1,2,isub)
                title('LAB Depth at 180 km')
            end
            
            ax=subplot(1,2,2);
            h3=ax.Children(5);
            h2=ax.Children(10);
            h1=ax.Children(15);
            leg=legend([h1,h2,h3],'10','20','40');
            leg.Location='none';
            leg.Box='off';
            leg.Title.String=({'LAB Amplitude' '(km)'});
            leg.Title.FontSize=7;
            
            leg.Units='Normalized';
            leg.Position=[0.6 0.65 0.1 0.1];
            

            
        end

    end
    methods (Access = private)
        function Inversion=ReadFromDisk(obj)
            Inversion=[];
            load(sprintf('%s/I',obj.path));
        end
    end    
end
