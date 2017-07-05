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
            
            k=kernel(3);
            do_column(2);
            newmodel=k.model();
            newmodel.hs=300;
            newmodel.vp=8.0;
            newmodel.vs=5.5;
            k=update_velocity_model(k,newmodel);
            do_column(1);
            
            function do_column(icol)
                ncol=2;
                
                isub=icol;
                ax1=subplot(4,2,isub);
                [G,XSP,A]=get_amplitude_components(k,0.10*111.11,-10);
                contourf(k.xs,k.zs,XSP');
                colorbar; hold on; colormap(ax1,parula)
                add_decor('Sp Scattering Pattern Contours');

                isub=isub+ncol;
                ax2=subplot(4,2,isub);
                contourf(k.xs,k.zs,log10(G)'); colorbar; hold on; colormap(ax2,parula);
                add_decor('log_1_0[ Geometrical Spreading Contours ]');

                isub=isub+ncol;
                ax3=subplot(4,2,isub);
                pcolor(k.xs,k.zs,A'); shading flat; colorbar; hold on; colormap(ax3,parula);
                add_decor('Timing Term');

                isub=isub+ncol;
                ax4=subplot(4,2,isub);
                K=XSP'.*G'.*A';
                pcolor(k.xs,k.zs,K); shading flat; colorbar; hold on; colormap(ax4,parula);
                add_decor('Product');
            end
            
        end
        
        function save(obj)
            obj.fig;
            print -depsc2 -painters figure4.eps
            close;
        end

    end
end
