classdef MyFigure
    properties
        fig
    end
    
    methods
        function obj=resize(obj,widthCm,heigthCm)
            set(obj.fig,'Units','centimeters')
            set(obj.fig,'Position',[0 0 widthCm heigthCm])
        end
        function save(obj,filename)
            obj.fig;
            tmpstr=sprintf('%s.eps',filename);
            print(tmpstr,'-depsc2', '-painters');
            close;
        end
    end
end