for depth = [120000 180000];
    for wl = [100000 200000 400000];
        for amp = [5000 10000 20000];
            savename=sprintf('%d-%d-%d',depth/1000,wl/1000,amp/1000);
            try
                FigureInversion(depth,wl,amp,savename);
            catch
                fprintf('Nothing found for %s\n',savename);
            end
        end
    end
end
