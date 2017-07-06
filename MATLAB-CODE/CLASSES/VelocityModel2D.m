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
        function plot_model(obj,varargin)
            clf;
            if nargin==2;
                filename=varargin{1};
            else
                filename='velocity_model';
            end

            figure(1);clf;
            set(gcf,'position',[0,0,800,800])
            subplot(2,3,[1 2 3]);

            pcolor(obj.xs,obj.zs,obj.dlnvs);
            shading flat
            title(sprintf('Conjugate-direction inversion\nVar Red = %.2f, nu=%.2e, norm opt = %d, Kernel Type: %d\n',obj.vred, obj.nu, obj.norm_opt, obj.Kernel_Type))
            c=colorbar;
            set(gca,'Ydir','reverse')
            ylabel('Depth [km]');
            xlabel('Offset [km]');
            %ylabel(c,'dlnv [\%]');
            ylabel(c,'(No Units Yet)')
            hold on;
            plot(obj.Locations,zeros(1,length(obj.Locations))+0.1,'o','markerfacecolor','black');
            plot(obj.xs,obj.xs*0+60,'--k')
            plot(obj.xs,120-5*cos(2*pi/400*(obj.xs-1725)),'--k');

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
                fprintf('iseis, isub = %d %d\n',iSeis,isub)

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
                    title(sprintf('%d of %d seismograms shown',iSeis,obj.nSeis));
                end
                xlabel('Time');
                axis off
            end

            tmpstr=sprintf('%s.eps',filename);
            print(tmpstr,'-depsc2', '-painters');
            close;
        end

        
    end
end