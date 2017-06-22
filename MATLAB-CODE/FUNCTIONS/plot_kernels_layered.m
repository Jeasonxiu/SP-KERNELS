function plot_kernels_layered()
%[Kernel,zs,xs,~,~] = analytical_kernel_layered();

Pdirect=0.09;
tchar=0.8;
nderiv=3.0;

model.hs=[60.0,60.0,180.0];
model.vp=[5.660,7.920,7.280];
model.vs=[3.2,4.4,4.05];

xs=(-600:1:150)';
zs=0:1:300;

%model.hs=[300];
%model.vp=[8.0];
%model.vs=[4.5];

[Kernel,~,~,~,~,~] = analytical_kernel_layered(Pdirect, xs, zs, tchar,nderiv, model);

[~,~,Pscat,~,~]=read_in_necessary_files(Pdirect, xs, zs, model);
[THE,~,~] = calc_angle(Pdirect*111.1,Pscat,xs,zs,model);

c=contourc(xs,zs,real(THE)',[0,1]);
s=getcontourlines(c);

npltx=2;
nplty=2;
nplt=npltx*nplty;

for iplt = 1:nplt;
   figure(1);
   subplot(npltx,nplty,iplt)
   pcolor(xs,zs,real(Kernel(:,:,iplt*25))); shading flat; colorbar; polarmap(); hold on
   plot3(s(1).x,s(1).y,zeros(1,length(s(1).x))+1,'-k'); hold on
   plot(0,0,'^r');
   ylabel('Depth (km)')
   xlabel('Offset (km)')
   set(gca,'YDir','reverse');
end

suptitle(sprintf('Sp Kernels in a Layered Medium, p = %.2f s/km',Pdirect))

%print('SpKern_Layers','-depsc2', '-painters');
%close;

end

function [THE,xs,zs] = calc_angle(Pdirect,Pscat,xs,zs,model)
    P1=Pdirect / 111.1;
    P2=Pscat / 111.1;
    
    %u
    [vps,vss]=get_v(model,zs);
    
    [XS,~]=meshgrid(xs,zs);
    [VPS,~]=meshgrid(vps,xs);
    [VSS,~]=meshgrid(vss,xs);
    
    THE1=asind(P1.*VSS);
    THE2=asind(P2.*VPS).*sign(-XS)';

    THE=THE1-THE2;
    
end

