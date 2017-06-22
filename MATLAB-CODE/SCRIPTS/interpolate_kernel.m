clear;
tmp=sprintf('Kernel_Angles_x2_%1.2f_%d',0.02,2500);
load([tmp '.mat']);

tDeci=1;

KTimes=KTimes(1:tDeci:end);

Kernel=Kernel(:,:,:,1:tDeci:end);

Scat_Depths_Fine=min(Scat_Depths):max(Scat_Depths);
Stations_Fine=1:0.2:150;
[STATIONS_FINE,SCAT_DEPTHS_FINE]=meshgrid(Stations_Fine,Scat_Depths_Fine);

Dims=size(Kernel);

Kernel_New=zeros(length(Scat_Depths_Fine),length(Stations_Fine),Dims(3),Dims(4));

for ii = 1:Dims(3)
    for jj = 1:Dims(4)
    Kernel_New(:,:,ii,jj)=interp2(Stations,Scat_Depths,Kernel(:,:,ii,jj),STATIONS_FINE,SCAT_DEPTHS_FINE,'cubic');
    end
end

Scat_Depths=Scat_Depths_Fine;
Stations=Stations_Fine;
Kernel=Kernel_New;

save([tmp '_interpolated.mat'],'-v7.3');