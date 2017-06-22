function plot_kernel_dependence_on_angle()

figure(1)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 300, 1000]); %<- Set size

load('Kernel_Angles_x2_0.02_2500_interpolated.mat')

stalocs_raw = load('stalocs.txt')/1000.0 - 1500.0;

stalocs = interp1(1:length(stalocs_raw),stalocs_raw,Stations);

itime=1;
iplot=0;

nangle=1;

dt=abs(KTimes(2)-KTimes(1));

for itime=20:20:100;
    for iangle = 1:nangle;
    
        figure(1)
        iplot=iplot+1;
        subplot(5,nangle,iplot);
        pcolor(-stalocs,-Scat_Depths,(Kernel(:,:,iangle,itime))); shading flat; caxis([-0.0002,0.0002]);
        ylim([-300,0]);
        title(sprintf('time, angle = %.0f s, %.0f^\\circ ',KTimes(itime),Angles(iangle)))
        hold on
        plot(0,-0.001,'marker','^','markerfacecolor','red')
        polarmap();
        add_isochron(KTimes(itime));
        if iplot == 21;
            ylabel('Scattering Depth [km]')
            xlabel('Scattering Offset [km]')
        end
        
        %plot the ray and fresnel zone
        c=4.5;
        T=4;
        wl=c*T;
        d=wl/3;
        deps=(0.0:300.0);
        xoffs=-deps/tand(90-Angles(iangle));
        plot(xoffs,-deps,'-black')
        FZHW=sqrt((10+deps).^2-deps.^2);
        plot(xoffs+FZHW,-deps,'--black')
        plot(xoffs-FZHW,-deps,'--black')
        
       
        
        
        ylabel('Depth (km)')
        colorbar;
   
    end
end

xlabel('Lateral Position (km)');

figure(1)
suptitle('Daughter kernels, positive weak scatterer')

set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 kernels_times.eps
close;

end

function add_isochron(time)
%
% Plots isochron on axiss
%
the=25.0;
vp=7.92;
vs=4.4;
[T,xs,zs]=timeshifts_halfspace(0,0,the,vp,vs);
C = contourc(xs,zs,T,[time,time]);
plot(C(1,2:end),-C(2,2:end),':k','LineWidth',1.5)


end

function [T,xs,zs] = timeshifts_halfspace(xrec,zrec,the,vp,vs)
%
%
%  vp=7.92;
%  vs=4.4;
%
% set up nodes
  x1=-600.0;
  x2=150.0;
  dx=5.0;
  
  z1=0.0;
  z2=300.0;
  dz=1.0;

  xs=x1:dx:x2;
  zs=z1:dz:z2;

  [XS,ZS]=meshgrid(xs,zs);
  
% calc traveltime for S-wave
  P1=[0.0 0.0];
  P2=[1.0 sind(the)];
  NUM= (P2(2)-P1(2))*XS - (P2(1)-P1(1))*ZS + P2(1)*P1(2) - P2(2)*P1(1);
  DEN=((P2(2)-P1(2)).^2 + (P2(1)-P1(1)).^2 )^0.5;
  
  T1=(NUM/DEN)/vs;
  
% calc time from node pt to scatterer

  T2=(sqrt((XS-xrec).^2+(ZS-zrec).^2))/vp;
  
  trec=interp2(XS,ZS,T1,xrec,zrec);
  
  T=T1+T2 - trec;

end
