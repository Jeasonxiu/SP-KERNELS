function [THE,xs,zs] = calc_angle(Pdirect,Pscat,xs,zs,model)
    P1=Pdirect / 111.1;
    P2=Pscat / 111.1;
    
    %u
    %[vps,vss]=get_velocity_from_profile('MIGRA/myvmod.nd',zs);
    
    [vps,vss]=get_v(model,zs);
    
    %[XS,~]=meshgrid(xs,zs);
    [VPS,~]=meshgrid(vps,xs);
    [VSS,~]=meshgrid(vss,xs);
    
    THE1=asind(P1.*VSS);
    THE2=asind(P2.*VPS);

    THE=THE1-THE2;
    
end