
path1='/Volumes/nmancine/data2/nmancine/PROJECTS/SP_RECEIVER_FUNCTIONS/KERNEL/SP-KERNELS/KERNEL-SEM/';

figure(1); clf; hold on;

iseis=75;

ipd=0;
for period = {'OUTPUT_FILES_200000-0.01-2500-15-DBLPERIOD', 'OUTPUT_FILES_200000-0.01-7500-15-DBLPERIOD'};
    ipd=ipd+1;
    path2=sprintf('%s/OUTPUT_FILES',period{1});
    isub=0;
    for channel = ['S','P'];
        isub=isub+1;
        subplot(2,1,isub);

        path3=sprintf('AA.S%04d.BX%s.semd',iseis,channel);

        path = sprintf('%s/%s/%s',path1,path2,path3);

        tr=load(path);

        if isub==2 && ipd==2;
            scl=(1/3)^2;
            style='--k';
        else
            scl=1.0*0.95;
            style='-r';
            
        end
        
        plot(tr(:,1),tr(:,2)*scl,style); hold on
        ylabel('u')
        xlabel('t')
        title(channel)
        %xlim([50 100])

    end

end

