function plot_model(filename,varargin)
    load(filename)
    clf;
    
    if nargin==2;
        ifig=varargin{1};
    else
        ifig=1;
    end
    
    figure(ifig);clf;
    set(gcf,'position',[0,0,800,800])
    subplot(2,3,[1 2 3]);

    pcolor(xs,zs,volume);
    shading flat
    title(sprintf('Conjugate-direction inversion\nVar Red = %.2f, nu=%.2e, norm opt = %d, Kernel Type: %d\n',vred,nu,norm_opt, Kernel_Type))
    c=colorbar;
    set(gca,'Ydir','reverse')
    ylabel('Depth [km]');
    xlabel('Offset [km]');
    %ylabel(c,'dlnv [\%]');
    ylabel(c,'(No Units Yet)')
    hold on;
    plot(Locations,zeros(1,length(Locations))+0.1,'o','markerfacecolor','black');
    plot(xs,xs*0+60,'--k')
    plot(xs,120-5*cos(2*pi/400*(xs-1725)),'--k');

    npts=length(dhat)/nSeis;  

    dt=KTimes(1)-KTimes(2);
    %time=(0:(npts-1)) * dt;
    time=KTimes;
    
    tmp1=1;
    tmp2=npts;
    
    nxplt=2;
    nyplt=6;
    nPerPlt=30;
    
    for iSeis = 1:1:nSeis;
        
        icount=0;
        isub=1;
        for ii = 1:iSeis;
            icount=icount+1;
            if (icount>nPerPlt)
                icount=1;
                isub=isub+1;
            end
        end
        fprintf('iseis, isub = %d %d\n',iSeis,isub)

        if isub>nyplt;
            iSeis=iSeis-1;
            break
        end
        
        subplot(nxplt,nyplt,nyplt+isub);

        scale=1./max(abs(d(tmp1:tmp2)))*1.5;

        plot(time,d(tmp1:tmp2)*scale+iSeis,'black')
        hold on
        plot(time,dhat(tmp1:tmp2)*scale+iSeis,'red')
        %title('Daughter Waveforms (every 10th one)')
        tmp1=tmp1+npts;
        tmp2=tmp2+npts;
        ylim([(isub-1)*nPerPlt-5,isub*nPerPlt+5])
        xlim([-20,0])

    end

    for jj = 1:nyplt;
        subplot(nxplt,nyplt,nyplt+jj)
        if jj==1;
            ylabel('Trace No.');
            handle=title(sprintf('%d of %d seismograms shown',iSeis,nSeis));
            %set(handle,'Position',[0.5,0.5]);
        end
        xlabel('Time');
        axis off
    end

    
    
    tmpstr=sprintf('%s.eps',filename);
    print(tmpstr,'-depsc2', '-painters');
    close;
end
