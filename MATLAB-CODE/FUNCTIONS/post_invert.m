   
function post_invert(filename)

    figure(1); clf;

    tmpstr=sprintf('%s.mat',filename);
    load(tmpstr)

    figure(1);subplot(2,3,[2 3 5 6]);
    pcolor(xs,zs,volume);
    title(sprintf('Nu = %e ',nu));
    shading flat
    title(sprintf(filename),'Interpreter','none')
    c=colorbar;
    set(gca,'Ydir','reverse')
    ylabel('Depth [km]');
    xlabel('Offset [km]');
    ylabel(c,'dlnv [%]');
    %caxis([-10 10]);
    hold on;
    plot(stalocs2,zeros(1,length(stalocs2))-0.01,'marker','^','markerfacecolor','red')
    
    nseis=150;
    
    npts=length(dhat)/nseis;
    
    dt=0.04*tDeci;
    time=(0:(npts-1)) * dt;
    
    dbest=B*m;
    
    tmp1=1;
    tmp2=npts;
    
    for iseis = 1:10:nseis;
    
        figure(1);subplot(2,3,[1 4]);

        scale=1./max(abs(dhat(tmp1:tmp2)))*5.0;

        plot(-time,dhat(tmp1:tmp2)*scale+iseis,'black')
        hold on
        plot(-time,dbest(tmp1:tmp2)*scale+iseis,'red')
        title('Daughter Waveforms (every 10th one)')
        tmp1=tmp1+npts;
        tmp2=tmp2+npts;
        
    end

    ylabel('Trace #');
    xlabel('Time');
    ylim([-2,102])
    
    tmpstr=sprintf('image_%s.png',filename);
    print(tmpstr,'-dpng','-r300');

end
