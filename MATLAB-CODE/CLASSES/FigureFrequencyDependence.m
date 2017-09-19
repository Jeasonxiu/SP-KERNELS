classdef FigureFrequencyDependence < MyFigure
	%Figure 3
    
    methods
        function obj=FigureFrequencyDependence()
            
            function add_decor()
                ylabel('Depth (km)');
                xlabel('Offset (km)');
                set(gca,'YDir','reverse')
                plot(0,0,'^r')
            end
            
            function add_subplot(isub,itime,ip,period)
                subplot(2,1,isub); hold on; box on;  
                pcolor(k.xs,k.zs,k.Kernel(:,:,ip,itime)); shading flat
                polarmap()
                add_decor();
                %title(sprintf('%.2f s, %2d degrees, %d-s period',k.KTimes(itime), k.Angles(ip), period),'Interpreter','Latex');
                title('')
                tmp=sprintf(...
                    '$t = %.1f$ s\n$\\theta  = %2d^\\circ$',k.KTimes(itime), k.Angles(ip));
                text(0.05,0.80,tmp,...
                    'Interpreter','Latex',...
                    'Units','Normalized',...
                    'HorizontalAlignment','left');
                title(sprintf('Period: %d s',period));
                
                
                xlim([-600,150]);
            end
            
            fig=figure(1);
            fig.Units='Inches';
            fig.Position=[0,0,4,8];
            
            k=kernel(1);
            k=load(k,'/Kernel_Angles_x2_0.01_2500_new.mat');
            obj.fig=figure(1);         
            add_subplot(1,80,3,4)
            xlabel('');
            k=kernel(1);
            k=load(k,'/Kernel_Angles_x2_0.01_2500.mat');
            add_subplot(2,80,3,8)
            
            
            
            
        end
        
        function save(obj)
            obj.fig;
            print -depsc2 -painters figureFrequencyDependence.eps
            close;
        end

    end
    
    methods (Static)
        function [Kernel,X,Y,Z,KTimes,nTimes,Angles] = load_kernel(path,tDeci)
            
            Kernel=[];
            KTimes=[];
            Angles=[];

            str1=[path '/Kernel_Angles_x2_0.01_7500_new.mat'];
            str2=[path '/stalocs.txt'];

            load(str1)
            stalocs = load(str2)/1000.0 - 1500.0; %convert to km

            %Decimate Kernels to speed up inversions
            %ip=1; %just use first angle calc for now

            %Going to time zero appears to work
            it1=1;

            Kernel=squeeze(Kernel(:,:,:,it1:tDeci:end));
            KTimes=KTimes(it1:tDeci:end);
            nTimes=length(KTimes);
            %Kernel=cumtrapz(Kernel,3);

            [X,Y,Z] = meshgrid(-stalocs,Scat_Depths,1:nTimes);
        end
        
    end
end
