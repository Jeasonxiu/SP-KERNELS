classdef figure3 < MyFigure
	%Figure 3
    
    methods
        function obj=figure3()
            
            function add_decor()
                ylabel('Depth (km)');
                xlabel('Offset (km)');
                set(gca,'YDir','reverse')
                plot(0,0,'^r')
            end
            
            function add_subplot(isub,itime,ip)
                subplot(2,3,isub); hold on; box on;  
                pcolor(k.xs,k.zs,k.Kernel(:,:,ip,itime)); shading flat
                polarmap()
                add_decor();
                title(sprintf('%.2f s, %2d deg',k.KTimes(itime), k.Angles(ip)));
                xlim([-600,150]);
            end
            
            k=kernel(1);
            k=load(k);
            obj.fig=figure(1);         
            add_subplot(1,120,1)
            add_subplot(2,80,1)
            add_subplot(3,40,1)
            add_subplot(4,80,1)
            add_subplot(5,80,2)
            add_subplot(6,80,3)
        end
        
        function save(obj)
            obj.fig;
            print -depsc2 -painters figure3.eps
            close;
        end

    end
end
