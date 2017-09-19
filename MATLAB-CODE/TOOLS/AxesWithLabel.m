classdef AxesWithLabel
    properties
        label
        Axes
    end
    methods
        function obj=AxesWithLabel(ax,label)
            obj.label=label;
            obj.Axes=ax;
            buf=0.05;
            tx=text(0,0,obj.label);
            tx.Units='Normalized';
            tx.FontSize=16;
            tx.FontWeight='Bold';
            tx.Position(1)=-4.5*buf;
            tx.Position(2)=1;
        end
    end
end