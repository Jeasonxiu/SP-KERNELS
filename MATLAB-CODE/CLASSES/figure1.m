classdef figure1
	%Figure 1
    properties
        fig
    end
    
    methods
        function obj=figure1()
            
            function add_decor()
                ylabel('Depth (km)');
                xlabel('Offset (km)');
                set(gca,'YDir','reverse')
                plot(0,0,'^r')
            end
            
            rp=0.09*111.11;
            %model.hs=[60.0,60.0,180.0];
            %model.vp=[5.660,7.920,7.280];
            %model.vs=[3.2,4.4,4.05];
            model.hs=[300];
            model.vp=[7.920];
            model.vs=[4.4];
            iso=isochrons(-600:15:150,0:2:300,rp,model);
            
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
