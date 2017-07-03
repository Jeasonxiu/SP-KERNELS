classdef figure3
	%Figure 3
    properties
        fig
    end
    
    methods
        function obj=figure3()
            
            function add_decor()
                ylabel('Depth (km)');
                xlabel('Offset (km)');
                set(gca,'YDir','reverse')
                plot(0,0,'^r')
            end
            
            kernel=kernel();
            
            obj.fig=figure(1); 
            subplot(2,1,1); hold on;
            [C,h]=contour(iso.xs,iso.zs,iso.Ps,0:6:40,'black');
            clabel(C,h,'FontWeight','bold','Color','blue')
            title('P-to-S');
            add_decor();
            
            subplot(2,1,2); hold on;
            [C,h]=contour(iso.xs,iso.zs,iso.Sp,-42:6:0,'black');
            clabel(C,h,'FontWeight','bold','Color','blue')
            title('S-to-P');
            add_decor();
        end
        
        function save(obj)
            obj.fig;
            print -depsc2 -painters figure1.eps
            close;
        end

    end
end
