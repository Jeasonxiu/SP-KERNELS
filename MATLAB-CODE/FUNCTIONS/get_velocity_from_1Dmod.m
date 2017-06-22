function [v]=get_velocity_from_1Dmod(model,z)
if z>sum(model.hs) || z<0
    fprintf('Error: Depth outside model range \n')
    v=-1;
    return
end
zref=0;
for ilay = 1:length(model.hs);
    ztop=zref;
    zbot=zref+model.hs(ilay);
    zref=zbot;
    if z <= zbot
        v=model.vs(ilay);
        return
    end
end

end

function test()

model.vs=[3.2,4.4,4.05];
model.hs=[60.0,60.0,180.0];

for z = -10:310;
    v=get_velocity_from_1Dmod(model,z);
    %fprintf('%f %f\n',z,v);
end

end
