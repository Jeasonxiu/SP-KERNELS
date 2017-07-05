classdef figure4
	%Figure 4
    properties
        fig
    end
    
    methods
        function obj=figure4()
            
            function add_decor(titlestring)
                ylabel('Depth (km)');
                xlabel('Offset (km)');
                set(gca,'YDir','reverse')
                title(titlestring);
                plot(0,0,'^r')
            end
            
            ax1=subplot(4,1,1);
            k=kernel(3);
            [G,XSP,A]=get_amplitude_components(k,0.10*111.11,-10);
            contourf(k.xs,k.zs,XSP');
            colorbar; hold on; colormap(ax1,parula)
            add_decor('Sp Scattering Pattern Contours');
            
            ax2=subplot(4,1,2);
            contourf(k.xs,k.zs,log10(G)'); colorbar; hold on; colormap(ax2,parula);
            add_decor('log_1_0[ Geometrical Spreading Contours ]');
            
            ax3=subplot(4,1,3);
            pcolor(k.xs,k.zs,A'); shading flat; colorbar; hold on; colormap(ax3,parula);
            add_decor('Timing Term');
            
            ax4=subplot(4,1,4);
            K=XSP'.*G'.*A';
            pcolor(k.xs,k.zs,K); shading flat; colorbar; hold on; colormap(ax4,parula);
            add_decor('Product');
            
        end
        
        function save(obj)
            obj.fig;
            print -depsc2 -painters figure4.eps
            close;
        end

    end
end
