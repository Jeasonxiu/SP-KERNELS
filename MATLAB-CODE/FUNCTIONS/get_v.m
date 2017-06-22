function [vp,vs]=get_v(model,zs)

    vp=zeros(1,length(zs))-999.0;
    vs=zeros(1,length(zs))-999.0;
    
    for iz = 1:length(zs);
        z=zs(iz);

        nlay=length(model.vp);
        zref=0;
        for ilay=1:nlay;
            top=zref;
            bot=zref+model.hs(ilay);
            zref=bot;

            if ( z>bot )    %depth below layer
               %
            elseif (z<top)  %depth above layer
               %
            else             %depth within layer
               vp(iz)=model.vp(ilay);
               vs(iz)=model.vs(ilay);
               break
            end

        end
    
    end

    return

end