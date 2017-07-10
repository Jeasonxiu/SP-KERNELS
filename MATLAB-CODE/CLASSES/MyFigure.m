classdef MyFigure
    properties
        fig
    end
    
    methods
        function obj=resize(obj,widthCm,heigthCm)
            set(obj.fig,'Units','centimeters')
            set(obj.fig,'Position',[0 0 widthCm heigthCm])
        end
    end
end