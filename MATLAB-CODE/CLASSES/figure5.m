classdef figure5
	%Figure 5
    properties
        fig
    end
    
    methods
        function obj=figure5()
            function add_isochron(xs,zs,time,incang)
                %
                % Plots isochron on axiss
                %
                the=incang;
                vp=7.92;
                vs=4.4;
                [T,xs,zs]=timeshifts_halfspace(xs,zs,0,0,the,vp,vs);
                C = contourc(xs,zs,T,[time,time]);
                plot(C(1,2:end),-C(2,2:end),':k','LineWidth',1.5)
            end
            function [T,xs,zs] = timeshifts_halfspace(xs,zs,xrec,zrec,the,vp,vs)

                [XS,ZS]=meshgrid(xs,zs);

                % calc traveltime for S-wave
                P1=[0.0 0.0];
                P2=[cosd(the) sind(the)];
                NUM= (P2(2)-P1(2))*XS - (P2(1)-P1(1))*ZS + P2(1)*P1(2) - P2(2)*P1(1);
                DEN=((P2(2)-P1(2)).^2 + (P2(1)-P1(1)).^2 )^0.5;

                T1=(NUM/DEN)/vs;

                % calc time from node pt to scatterer

                T2=(sqrt((XS-xrec).^2+(ZS-zrec).^2))/vp;

                trec=interp2(XS,ZS,T1,xrec,zrec);

                T=T1+T2 - trec;

            end

            AddColorbars=false;
            
            KSEM_NFO=kernel(1);
            KSEM_NFO=load(KSEM_NFO);
            
            KAN_NFO=kernel(3);
            KAN_NFO.xs=KSEM_NFO.xs;
            KAN_NFO.zs=KSEM_NFO.zs;
            KAN_NFO.Angles=KSEM_NFO.Angles;
            newmodel=KAN_NFO.model;
            newmodel.vp=KAN_NFO.model.vp(2);
            newmodel.vs=KAN_NFO.model.vs(2);
            newmodel.hs=300;
            KAN_NFO=update_velocity_model(KAN_NFO,newmodel);
            KAN_NFO=load(KAN_NFO);

            it=2;
            itimes=[30,80,110];

            obj.fig=figure(1); clf;
            set(obj.fig,'defaulttextinterpreter','latex');
            set(obj.fig,'PaperPositionMode','auto')
            set(obj.fig,'DefaultAxesFontSize',6);
            xlimits=[-600, 150];
            for iangle=1:length(KSEM_NFO.Angles);
                itime=itimes(it);
                KSEM=flatten(KSEM_NFO.Kernel(:,:,iangle,itime),2);

                %Normalize by area
                KSEM=KSEM/(pi*2.5^2);

                iplt=iangle;

                subplot(3,3,0+iplt);
                pcolor(KSEM_NFO.xs,-KSEM_NFO.zs,KSEM); shading flat; hold on
                text(0.05,0.80,sprintf('SEM\nTime: %.2f s\nAngle: %.2f$^\\circ$', KSEM_NFO.KTimes(itime), KSEM_NFO.Angles(iangle)),'Units','normalized','Interpreter','Latex');
                polarmap();
                if AddColorbars;
                    colorbar;
                    clim=max(max(abs(KSEM)));
                    caxis([-clim,clim]);
                end

                xlim(xlimits)
                add_isochron(KSEM_NFO.xs,KSEM_NFO.zs,KSEM_NFO.KTimes(itime),KSEM_NFO.Angles(iangle))

                AnalyticalVersion=1;
                if AnalyticalVersion==-1;

                    %Angles=[15.0,20.0,25.0];
                    %[Kernel2,zs,xs,Angles,KTimes] = analytical_kernel_halfspace(Angles);
                    %KAN=flatten(Kernel2(:,:,iangle,itime),2);
                else
                    KAN=flatten(KAN_NFO.Kernel(:,:,itime,iangle),2);
                end

                %Normalize by area
                %fprintf('Normalizing by area\n')
                deltax=abs(KAN_NFO.xs(2)-KAN_NFO.xs(1));
                deltaz=abs(KAN_NFO.zs(2)-KAN_NFO.zs(1));
                KAN=KAN/(deltax*deltaz);

                subplot(3,3,3+iplt);
                pcolor(KAN_NFO.xs,-KAN_NFO.zs,KAN); shading flat; hold on;

                if AddColorbars;
                    colorbar;
                    clim=max(max(abs(KAN)));
                    caxis([-clim,clim]);
                end


                text(0.05,0.80,sprintf('Ray Theory\nTime: %.2f s\nAngle: %.2f$^\\circ$', KAN_NFO.KTimes(itime), KAN_NFO.Angles(iangle)),'Units','normalized','Interpreter','Latex');
                %xlabel('Lateral Position (km)')
                %ylabel('Depth (km)')
                xlim(xlimits)
                add_isochron(KAN_NFO.xs,KAN_NFO.zs,KAN_NFO.KTimes(itime),KAN_NFO.Angles(iangle))

                subplot(3,3,6+iplt)
                %pcolor(xs,-zs,fliplr(KSEM) - KAN/173.4297); shading flat; colorbar;
                %title(sprintf('Time, Angle = %f,  %f', KTimes(itime), Angles(iangle)))
                tmpf1=mean(abs(KAN));
                plot(KAN_NFO.xs,tmpf1); hold on;
                tmpf2=mean(abs(KSEM));
                scalingfac=tmpf1'\fliplr(tmpf2)';

                plot(KAN_NFO.xs,tmpf2/scalingfac,'-r')
                xlim(xlimits)
                text(0.05,0.85,sprintf('SEM/Analytic = %f',scalingfac),'Units','normalized');

            end

            subplot(3,3,1)
            ylabel('Depth (km)')
            subplot(3,3,4)
            ylabel('Depth (km)')
            subplot(3,3,7)
            ylabel('Mean Abs. Amplitude')
            xlabel('Lateral Position (km)')
            subplot(3,3,8)
            xlabel('Lateral Position (km)')
            subplot(3,3,9)
            xlabel('Lateral Position (km)')

            %FIG.PaperUnits='inches';
            %FIG.PaperPosition=[0 0 20 20];

            set(findall(gcf,'type','text'),'FontSize',6);
        
        end

        function save(obj)
            obj.fig;
            print -depsc2 -painters figure5.eps
            close;
        end

    end
end
