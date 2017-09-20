classdef VelocityModel2D
    properties
        dlnvs
        dlnvp
        dlnrho
        xs
        zs
        d
        dhat
        dtime
        vred
        nu
        norm_opt
        Kernel_Type
        nSeis
        Locations
    end
    methods
        function obj=VelocityModel2D()
        end
        function fig=plot_model(obj,varargin)
            set(0,'defaulttextInterpreter','latex') %latex axis labels
            
            if nargin>=2;
                filename=varargin{1};
            else
                filename='velocity_model';
            end
            
            if nargin>=3;
                clabel=varargin{2};
            else
                clabel='No Units';
            end
            
            if nargin>=4;
                labwavlen=varargin{3}.wavlen;
                labamp=varargin{3}.amplitude;
                labdep=varargin{3}.depth;
            end
            
            if nargin>=5;
                fig=gcf;
                %set(gcf,'position',[0,0,800,800])
                ax1=varargin{4};
                ImageOnly=true;
            else
                fig=figure(1);clf;
                %set(gcf,'position',[0,0,800,800])
                ax1=subplot(2,6,[2 3 4 5 6]);
                ImageOnly=false;
            end
            
            %get background profile
            %VM=velocity_model;
            %profile=get_v(VM,obj.zs);
            %[bkgrnd,~]=meshgrid(profile,obj.xs);
            %tmp=bkgrnd'.*(1+obj.dlnvs/100.0);

            pcolor(obj.xs,obj.zs,obj.dlnvs);
            shading flat
            polarmap()
            %title(sprintf('Conjugate-direction inversion\nVar Red = %.2f, nu=%.2e, norm opt = %d, Kernel Type: %d\n',obj.vred, obj.nu, obj.norm_opt, obj.Kernel_Type))

            
            set(gca,'Ydir','reverse')
            ylabel('Depth [km]');
            xlabel('Offset [km]');
            %ylabel(c,'dlnv [\%]');
            
            hold on;
            %plot(obj.Locations,zeros(1,length(obj.Locations))+0.1,'od','markerfacecolor','black','markersize',2,'markeredgecolor','none');
            l=plot(obj.Locations,zeros(1,length(obj.Locations)),'.','markersize',3);
            l.MarkerSize=1;
            l.MarkerEdgeColor='black';
            l.MarkerFaceColor='black';

            
            
            deltasta=obj.Locations(2)-obj.Locations(1);
            titstr=sprintf('Station spacing: $\\sim %d$ km\n',round(deltasta));
            title(titstr,'interpreter','latex');
            
            if nargin >= 4;
                plot(obj.xs,obj.xs*0+60,'--k')

                wl=labwavlen/1000.0;
                amp=labamp/1000.0;
                dep=labdep/1000.0;

                plot(obj.xs,dep-amp*cos(2*pi/wl*(obj.xs-1725)),'--k');
                
                if (ImageOnly==false)
                    titstr=[ sprintf('LAB Properties\n')...
                             sprintf('Amplitude (Trough-to-Peak) = %d km\n', amp*2.0)...
                             sprintf('$\\lambda/2$ = %d km\n', wl/2.0)...
                             sprintf('Depth = %d km', dep)];
                    buf=0.0125;

                    text(buf,1-buf,titstr,...
                        'interpreter','latex',...
                        'Units','normalized',...
                        'horizontalAlignment', 'left',...
                        'verticalAlignment', 'top');
                end
                
            end
            
            if ImageOnly;
                return
            end
            
            ax2=subplot(2,6,1); cla; 
            %plot(mean(obj.dlnvs'),obj.zs); hold on
            %plot(median(obj.dlnvs'),obj.zs)
            imid1=ceil(length(obj.xs)*1/3);
            imid2=ceil(length(obj.xs)/2);
            imid3=ceil(length(obj.xs)*2/3);
            plot(obj.dlnvs(:,imid1),obj.zs,'-'); hold on
            plot(obj.dlnvs(:,imid2),obj.zs,'-.')
            plot(obj.dlnvs(:,imid3),obj.zs,'--')
            set(gca,'Ydir','reverse')
            ylabel('Depth [km]')
            
            function minmax(imid,zmin,zmax)
                scatpot=obj.dlnvs(:,imid);
                
                vmax=min(scatpot);
                vmin=max(scatpot);
                
                for iz = 1:length(scatpot);
                    if obj.zs(iz) >= zmin && obj.zs(iz) <= zmax
                        vmax=max([vmax scatpot(iz)]);
                        vmin=min([vmin scatpot(iz)]);
                    end
                end
                
                fprintf('imid, vmin, vmax, diff = %d, %f, %f, %f\n',...
                    imid, vmin, vmax, vmax-vmin)
                
            end
            
            
            for imid = [imid1 imid2 imid3];
                minmax(imid,0,100);
                minmax(imid,100,250);
            end
            

            ax1.YLabel.String='';
            ax2.YTick=ax1.YTick;
            ax1.YTick=[];
            %set(ax1,'YLabel','none')
            
            

            npts=length(obj.dhat)/obj.nSeis;
            
            tmp1=1;
            tmp2=npts;

            nxplt=2;
            nyplt=6;
            nPerPlt=30;

            for iSeis = 1:1:obj.nSeis;

                icount=0;
                isub=1;
                for ii = 1:iSeis;
                    icount=icount+1;
                    if (icount>nPerPlt)
                        icount=1;
                        isub=isub+1;
                    end
                end
                %fprintf('iseis, isub = %d %d\n',iSeis,isub)

                if isub>nyplt;
                    iSeis=iSeis-1;
                    break
                end

                subplot(nxplt,nyplt,nyplt+isub);
                
                scale=1./max(abs(obj.d(tmp1:tmp2)))*1.5;

                plot(obj.dtime,obj.d(tmp1:tmp2)*scale+iSeis,'black')
                hold on
                plot(obj.dtime,obj.dhat(tmp1:tmp2)*scale+iSeis,'red')
                %title('Daughter Waveforms (every 10th one)')
                tmp1=tmp1+npts;
                tmp2=tmp2+npts;
                ylim([(isub-1)*nPerPlt-5,isub*nPerPlt+5])
                xlim([-20,0])

            end

            for jj = 1:nyplt;
                subplot(nxplt,nyplt,nyplt+jj)
                if jj==1;
                    ylabel('Trace No.');
                    ax2=gca();
                    axes(ax1);
                    ax1=text(0.5,-0.4,...
                        sprintf('%d of %d seismograms shown',iSeis,obj.nSeis),...
                        'interpreter','latex','units','normalized',...
                        'horizontalAlignment','center');
                    axes(ax2);
                    
                end
                xlabel('Time');
                axis off
            end

            tmpstr=sprintf('%s.eps',filename);
            print(tmpstr,'-depsc2', '-painters');
            close;
        end
        
        function [therat,medrat,stdrat,quant10,quant90]=plot_recovery_performance(obj)
            
            function [vswing] = minmax(imid,zmin,zmax)
                scatpot=obj.dlnvs(:,imid);
                
                vmax=min(scatpot);
                vmin=max(scatpot);
                
                for iz = 1:length(scatpot);
                    if obj.zs(iz) >= zmin && obj.zs(iz) <= zmax
                        vmax=max([vmax scatpot(iz)]);
                        vmin=min([vmin scatpot(iz)]);
                    end
                end
                
                vswing=vmax-vmin;
                %fprintf('imid, vmin, vmax, diff = %d, %f, %f, %f\n',...
                %    imid, vmin, vmax, vmax-vmin)
                
            end
            
            function [vrms] = getvrms(imid,zmin,zmax)
                scatpot=obj.dlnvs(:,imid);
                
                vsumsq=0;
                n=0;
                for iz = 1:length(scatpot);
                    if obj.zs(iz) >= zmin && obj.zs(iz) <= zmax
                        vsumsq=vsumsq+scatpot(iz)^2;
                        n=n+1;
                    end
                end
                vmeansq=vsumsq/n;
                vrms=sqrt(vmeansq);
        
                
                %fprintf('imid, vmin, vmax, diff = %d, %f, %f, %f\n',...
                %    imid, vmin, vmax, vmax-vmin)
                
            end
            
            %fig=figure(1);
            %fig.Units='Inches';
            %fig.Position(1)=0;
            %fig.Position(2)=0;
            %fig.Position(3)=6;
            %fig.Position(4)=8;
            
            
            for iter = 1:2;
                if iter == 1
                    offset_min=0;
                    offset_max=999999;
                    ls='--';
                else
                    offset_min=1320;
                    offset_max=2125;
                    ls='-';
                end

                icount=0;
                
                offsets= obj.xs( obj.xs >= offset_min & obj.xs <= offset_max );
                indeces= 1:length(obj.xs);
                indeces= indeces( obj.xs >= offset_min & obj.xs <= offset_max );
                moho_brightness=offsets*0;
                lab_brightness=offsets*0;
                
                for icol = 1:length(offsets);
                    icount=icount+1;
                    vswing=minmax(indeces(icol),10,100);
                    %vrms=getvrms(indeces(icol),10,100);
                    moho_brightness(icount)=vswing;
                    vswing=minmax(indeces(icol),110,200);
                    %vrms=getvrms(indeces(icol),110,200);
                    lab_brightness(icount)=vswing;
                end

                %subplot(5,1,1)
                %plot(offsets,moho_brightness,'LineWidth',1,'color','k','Linestyle',ls); hold on
                %ylabel({'RMS','Moho Amplitude'})
                %subplot(5,1,2)
                %plot(offsets,lab_brightness,'LineWidth',1,'color','k','Linestyle',ls); hold on
                %ylabel({'RMS','LAB Amplitude'})

                %xlabel('Lateral Offset (km)')
            end
            
            
            
            %subplot(9,1,[5:9])
            
            ratio=lab_brightness./median(moho_brightness);
            

           
            
            %plot(lab_brightness,moho_brightness,'+r'); hold on
            %xlabel('RMS LAB Amplitude')
            %ylabel('RMS Moho Amplitude')
            
            vref=velocity_model().vs(2);
            dlnv1=(velocity_model().vs(1)/velocity_model().vs(2))-1;
            dlnv3=(velocity_model().vs(2)/velocity_model().vs(3))-1;
            
            coef = abs(dlnv1/dlnv3);
            
            medrat=median(ratio);
            stdrat=std(ratio);
            quant10=quantile(ratio,0.10);
            quant90=quantile(ratio,0.90);
            therat=1./coef;
            
            return
            
            
            
            xtmp=15:40;
            ytmp=100:140;
            [XTMP,YTMP]=meshgrid(xtmp,ytmp);
            err = YTMP./(coef*XTMP) - 1;
            contours=-0.5:0.25:1.75;
            %[C,h]=contour(xtmp,ytmp,err,contours,'color','black');
            %clabel(C,h,'Fontsize',15,'Color','blue','LabelSpacing',500);
            
            [C,h]=contour(xtmp,ytmp,err,contours,'color','black');
            tl=clabel(C,h,'manual','EdgeColor','white', 'BackgroundColor','white');
            
            for ii =1:length(tl)
                newstr=sprintf('%d %%',str2num(tl(ii).String)*100);
                tl(ii).String=newstr;
            end
            
        end
             
    end
end