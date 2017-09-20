classdef figure3 < MyFigure
	%Figure 3
    
    methods
        function obj=figure3()
            

            
            function add_subplot(isub,itime,ip)
                
                function add_decor()
                    ylabel('Depth (km)');
                    xlabel('Offset (km)');
                    set(gca,'YDir','reverse')
                    plot(0,0,'^r')
                    c=colorbar();

                    ax.Position(4)=ax.Position(4)*1.0;
                    totwid=ax.Position(3);

                    gr=1.618;
                    cwid=totwid*0.03;
                    buf=cwid*gr;

                    ax.Position(3)=totwid*1.00;
                    c.Position(1)=ax.Position(1)+totwid+buf;
                    c.Position(2)=ax.Position(2);
                    c.Position(3)=cwid;
                    c.Position(4)=ax.Position(4)/gr;
             
                
                end
                
                
                ax=subplot(2,3,isub); hold on; box on;  
                pcolor(k.xs,k.zs,k.Kernel(:,:,ip,itime)); shading flat
                polarmap()
                add_decor();
                tmp=sprintf(...
                    '$t = %.1f$ s\n$\\theta  = %2d^\\circ$',k.KTimes(itime), k.Angles(ip));
                text(0.05,0.80,tmp,...
                    'Interpreter','Latex',...
                    'Units','Normalized',...
                    'HorizontalAlignment','left');
                xlim([-600,150]);
                
                if isub ~= 1 && isub ~= 4
                    ax.YTickLabel={};
                    ylabel('')
                end
                if isub <= 3
                    xlabel('')
                end
                
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
