function [T,X]=shoot_ray(p,z,model)
vs=model.vs;
hs=model.hs;

nlay=length(hs);
T=0.0;
X=0.0;
zref=0.0;
for ilay=1:nlay;
    top=zref;
    bot=zref+model.hs(ilay);

    if ( z>=bot )    %depth below layer
        dz=0.0;
    elseif (z<=top)  %depth above layer
        dz=hs(ilay);
    else             %depth within layer
        dz=hs(ilay)-(z-top);    
    end
        
    zref=bot;
    u=1./vs(ilay);
    T=T+u.^2*dz/(u.^2-p.^2).^0.5;
    X=X+p.*dz/(u.^2-p.^2).^0.5;
    
end

end

function calculate_travel_time_field()

xs=-600:150;
zs=0:300;

model.vs=[3.2,4.4,4.05];
model.hs=[60.0,60.0,180.0];

p=0.09;

Tss=zeros(1,length(zs));
Xss=zeros(1,length(zs));
i=0;
for z = zs;
    i=i+1;
    [T,X]=shoot_ray(p,z,model);
    Tss(i)=T;
    Xss(i)=X;
    Tss(i)=T-p*X;
    fprintf('%f %f %f %f\n',z,T,X,Tss(i))
end
%clf;
%plot(Xss,Tss);

[XS,TSS]=meshgrid(xs,Tss);
TT = TSS + XS*p;

contourf(xs,-zs,TT); colorbar;

end