clear;

figure(1)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1200, 1000]); %<- Set size

figure(2)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1200, 1000]); %<- Set size

figure(3)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1200, 1000]); %<- Set size

load Kernel_Plus_Weaker

tmp1=load('Kernel_Plus_Weaker','Kernel');
tmp2=load('Kernel_Minus_Weaker','Kernel');

Kernel1=tmp1.Kernel;
Kernel2=tmp2.Kernel;
KernelSum=Kernel1+Kernel2;

stalocs = load('stalocs.txt')/1000.0;

for iplot = 1:1:25;
    
    itime=iplot*30;
    
    figure(1)
    subplot(5,5,25+1-iplot);
    pcolor(-stalocs+1650,-Scat_Depths,(Kernel1(:,:,itime))); caxis([-0.001,0.001]); colorbar; shading flat;
    title(sprintf('RF Amp. at Time = %.0f ',KTimes(itime)))
    hold on
    plot(0,-10.001,'marker','^','markerfacecolor','red')
    polarmap();
    %
    figure(2)
    subplot(5,5,25+1-iplot);
    pcolor(-stalocs+1650,-Scat_Depths,(Kernel2(:,:,itime))); caxis([-0.001,0.001]); colorbar; shading flat;
    title(sprintf('RF Amp. at Time = %.0f ',KTimes(itime)))
    hold on
    plot(0,-10.001,'marker','^','markerfacecolor','red')
    polarmap();
    %
    figure(3)
    subplot(5,5,25+1-iplot);
    pcolor(-stalocs+1650,-Scat_Depths,(KernelSum(:,:,itime))); caxis([-0.001,0.001]); colorbar; shading flat;
    title(sprintf('RF Amp. at Time = %.0f ',KTimes(itime)))
    hold on
    plot(0,-10.001,'marker','^','markerfacecolor','red')
    polarmap();
end

figure(1)
suptitle('Positive scatterer')
print('kernels_positive','-dpng','-r300');

figure(2)
suptitle('Negative scatterer')
print('kernels_negative','-dpng','-r300');

figure(3)
suptitle('Positive + Negative')
print('kernels_sum','-dpng','-r300');



