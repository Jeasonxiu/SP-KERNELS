classdef figure1
	%figure 1
    
    methods
        function obj=figure1()
            
            function add_decor()
                ylabel('Depth (km)');
                xlabel('Offset (km)');
                set(gca,'YDir','reverse')
                plot(0,0,'^r')
            end
            
            rp=0.09*111.11;
            
            iso=isochrons(-600:15:150,0:2:300,rp);
            
            figure(1); 
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

    end
end
