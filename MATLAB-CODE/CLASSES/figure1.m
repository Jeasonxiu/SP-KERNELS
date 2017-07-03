classdef figure1
	%figure 1
    properties
        isochrons=zeros(50,50);
    end
    
    methods 
        function plot(obj)
            function add_labels()
                ylabel('Depth (km)');
                xlabel('Offset (km)');
            end
            
            figure(1); hold on;
            subplot(2,1,1)
            pcolor(obj.isochrons);
            title('P-to-S');
            add_labels();
            
            subplot(2,1,2)
            pcolor(obj.isochrons);
            title('S-to-P');
            add_labels();
            
        end

    end
end
