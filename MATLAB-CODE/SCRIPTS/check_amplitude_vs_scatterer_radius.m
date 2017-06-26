%A test to see how a 1% pct scatterer changes amplitude with size.
%
% Normalization by the scatterer area seems to work for radii > 5 km.
%
%
%%loop over radii to compare
figure(1);clf;
radii=[2500, 5000, 7500];
i=0;
for radius = radii ;
    seis=load(sprintf('../../../../OUTPUT_FILES_95000-0.01-%d-15-DBLPERIOD/OUTPUT_FILES/AA.S0020.BXP.semd',radius));
    
    area=pi*radius^2;
    fac=1.0/area;
    
    plot(seis(:,1),seis(:,2)*fac); hold on;
    
    i=i+1;
    labels{i}=num2str(radius);
    
end
legend(labels)
    
    