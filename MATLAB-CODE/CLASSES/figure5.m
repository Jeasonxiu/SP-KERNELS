classdef figure5
	%Figure 5
    properties
        fig
    end
    
    methods
        function obj=figure5()
            
        end
        
        function save(obj)
            obj.fig;
            print -depsc2 -painters figure5.eps
            close;
        end

    end
end
