function [m]=mask_out_model(m,nOff,nDep,nSeis,G,dhat)

    icol=0;
    for ioff = 1:1:nOff;
        for idep = 1:1:nDep;
           icol=icol+1;
           volume(idep,ioff)=m(icol);
           if idep<50; volume(idep,ioff)=0.0; end;
        end
    end
    
    pcolor(volume); shading flat;
    
    icol=0;
    for ioff = 1:1:nOff;
        for idep = 1:1:nDep;
           icol=icol+1;
           m(icol)=volume(idep,ioff); 
        end
    end
    
    
    npts=length(dhat)/nSeis;

    dt=KTimes(1)-KTimes(2);
    %time=(0:(npts-1)) * dt;
    time=KTimes;
    
    tmp1=1;
    tmp2=npts;
    for iseis = 1:1:nSeis;

        if iseis<=20;
            isub=0;
        elseif iseis<=40;
            isub=1;
        elseif iseis<=60;
            isub=2;
        else
            continue
        end

        subplot(2,3,4+isub);

        scale=1./max(abs(d(tmp1:tmp2)))*1.5;

        plot(time,d(tmp1:tmp2)*scale+iseis,'black')
        hold on
        plot(time,dhat(tmp1:tmp2)*scale+iseis,'red')
        %title('Daughter Waveforms (every 10th one)')
        tmp1=tmp1+npts;
        tmp2=tmp2+npts;

    end

    subplot(2,3,4)
    ylabel('Trace #');
    xlabel('Time');
    %ylim([-2,102]);
    
    
end