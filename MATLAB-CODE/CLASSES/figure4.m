classdef figure4 < MyFigure
	%Figure 4
    
    methods
        function obj=figure4()
            

            
            obj.fig=figure(1);
            set(0,'defaulttextinterpreter','latex')
            k=kernel(3);
            do_column(2);
            newmodel=k.model;
            newmodel.hs=300;
            newmodel.vp=k.model.vp(2);
            newmodel.vs=k.model.vs(2);
            k=update_velocity_model(k,newmodel);
            do_column(1);
            
            function do_column(icol)
                function add_decor(titlestring)
                    ylabel('Depth (km)');
                    xlabel('Offset (km)');

                    set(gca,'YDir','reverse')
                    set(gca,'Fontsize',7)
                    title(titlestring);
                    plot(0,0,'^r')

                    if icol==2; ylabel(''); end
                    if isub<=6; xlabel(''); end
                    
                    
                    ax.Position(4)=ax.Position(4)*0.95;
                    totwid=ax.Position(3);
                    
                    gr=1.618;
                    cwid=totwid*0.03;
                    buf=cwid*gr;
                    
                    ax.Position(3)=totwid*1.00;
                    c.Position(1)=ax.Position(1)+totwid+buf;
                    c.Position(2)=ax.Position(2);
                    c.Position(3)=cwid;
                    c.Position(4)=ax.Position(4)/gr;
                    
                    labs={'(a)','(e)','(b)','(f)','(c)','(g)','(d)','(h)'};
                    lab=labs{isub};
                    AxesWithLabel(ax,lab);

                end             
                
                ncol=2;
                
                polar=polarmap();
                
                isub=icol;
                ax=subplot(4,2,isub);
                [G,XSP,A]=get_amplitude_components(k,0.10*111.11,-10);
                contourf(k.xs,k.zs,XSP');
                c=colorbar; hold on; colormap(ax,polar);
                caxis([-1,1]*max(abs(caxis)));
                add_decor('Sp Scattering Pattern Contours');
                
                isub=isub+ncol;
                ax=subplot(4,2,isub);
                contourf(k.xs,k.zs,log10(G)'); c=colorbar; hold on;
                colormap(ax,parula);
                add_decor('$\log_{10}$[ Geometrical Spreading Contours ]');
                
                isub=isub+ncol;
                ax=subplot(4,2,isub);
                pcolor(k.xs,k.zs,A'); shading flat; c=colorbar; hold on;
                colormap(ax,polar); caxis([-1,1]*max(abs(caxis)));
                add_decor('Timing--Wavelet Term');

                isub=isub+ncol;
                ax=subplot(4,2,isub);
                K=XSP'.*G'.*A';
                pcolor(k.xs,k.zs,K); shading flat; c=colorbar; hold on;
                colormap(ax,polar); caxis([-1,1]*max(abs(caxis)));
                add_decor('Product');
                polarmap;
                

                
                
                
            end
            
        end
        
        function save(obj)
            obj.fig;
            print -depsc -painters figure4.eps
            close;
        end

    end
end
