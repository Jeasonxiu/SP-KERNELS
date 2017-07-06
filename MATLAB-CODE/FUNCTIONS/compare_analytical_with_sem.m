function compare_analytical_with_sem()
    clear;
    
    AddColorbars=false;

    load('/Volumes/nmancine/data2/nmancine/PROJECTS/SP_RECEIVER_FUNCTIONS/KERNEL/SP-KERNELS/DATA/Kernel_Angles_x2_0.01_2500.mat');

    stalocs_raw = load('/Volumes/nmancine/data2/nmancine/PROJECTS/SP_RECEIVER_FUNCTIONS/KERNEL/SP-KERNELS/DATA/stalocs.txt')/1000.0 - 1500.0;

    stalocs = interp1(1:length(stalocs_raw),stalocs_raw,Stations);

    clearvars -except Kernel KTimes stalocs Scat_Depths Angles AddColorbars

    it=2;
    itimes=[30,80,110];

    FIG=figure(1); clf;
    set(FIG,'defaulttextinterpreter','latex');
    set(FIG,'PaperPositionMode','auto')
    set(FIG,'DefaultAxesFontSize',6);
    xlimits=[-600, 150];
    for iangle=1:3;
        itime=itimes(it);
        KSEM=flatten(Kernel(:,:,iangle,itime),2);
        
        %Normalize by area
        KSEM=KSEM/(pi*2.5^2);

        iplt=iangle;

        subplot(3,3,0+iplt);
        pcolor(-stalocs,-Scat_Depths,KSEM); shading flat; hold on
        text(0.05,0.80,sprintf('SEM\nTime: %.2f s\nAngle: %.2f$^\\circ$', KTimes(itime), Angles(iangle)),'Units','normalized','Interpreter','Latex');
        polarmap();
        if AddColorbars;
            colorbar;
            clim=max(max(abs(KSEM)));
            caxis([-clim,clim]);
        end
        
        xlim(xlimits)
        add_isochron(-stalocs,Scat_Depths,KTimes(itime),Angles(iangle))

        AnalyticalVersion=1;
        if AnalyticalVersion==-1;
            
            Angles=[15.0,20.0,25.0];
            [Kernel2,zs,xs,Angles,KTimes] = analytical_kernel_halfspace(Angles);
            KAN=flatten(Kernel2(:,:,iangle,itime),2);
        else
        %%%%%New section
            %model.hs=300.0;
            %model.vp=7.92;
            %model.vs=4.4;

            model=velocity_model();
            newmodel=model;
            newmodel.hs=300.0;
            newmodel.vp=7.92;
            newmodel.vs=4.4;
            model=update(model,newmodel);

            Angles=[15.0,20.0,25.0];
            Pdirect=sind(Angles)/model.vs;
            tchar=1.6;
            nderiv=3.0;

            xs=linspace(-600,150,150)';
            zs=0:5:300;

            [Kernel2,~,~,~,~,~] = analytical_kernel_layered(Pdirect, xs, zs, tchar,nderiv, model);
            KAN=flatten(Kernel2(:,:,itime,iangle),2);
        end
        
        %Normalize by area
        %fprintf('Normalizing by area\n')
        deltax=xs(2)-xs(1);
        deltaz=zs(2)-zs(1);
        KAN=KAN/(deltax*deltaz);
        
        %debug
        fprintf('deltax, deltaz = %10f, %10f', deltax, deltaz )

        subplot(3,3,3+iplt);
        pcolor(xs,-zs,KAN); shading flat; hold on;
        
        if AddColorbars;
            colorbar;
            clim=max(max(abs(KAN)));
            caxis([-clim,clim]);
        end

        
        text(0.05,0.80,sprintf('Ray Theory\nTime: %.2f s\nAngle: %.2f$^\\circ$', KTimes(itime), Angles(iangle)),'Units','normalized','Interpreter','Latex');
        %xlabel('Lateral Position (km)')
        %ylabel('Depth (km)')
        xlim(xlimits)
        add_isochron(xs,zs,KTimes(itime), Angles(iangle))

        subplot(3,3,6+iplt)
        %pcolor(xs,-zs,fliplr(KSEM) - KAN/173.4297); shading flat; colorbar;
        %title(sprintf('Time, Angle = %f,  %f', KTimes(itime), Angles(iangle)))
        tmpf1=mean(abs(KAN));
        plot(xs,tmpf1); hold on;
        tmpf2=mean(abs(KSEM));
        scalingfac=tmpf1'\fliplr(tmpf2)';

        plot(-stalocs,tmpf2/scalingfac,'-r')
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

    !rm kernel_comparison.eps
    print -depsc2 -painters kernel_comparison.eps
    close;

    end

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
