classdef ReceiverFunction
    properties
        Parent
        Daughter
        rf
        time
    end
    
    methods
        function obj=ReceiverFunction(Parent,Daughter)
 
            
            
            [~,imax]=max(abs(Parent.u));
            
            Mask=(Parent.t(:)>Parent.t(imax-40)).*(Parent.t(:)<Parent.t(imax+40));
            Parent.u=Parent.u.*Mask;
            
            TB=4;
            NT=3;
            dt=Parent.dt;
            win_len=100;
            Poverlap=0.90;

            it_end=imax+round(win_len/2/dt);
            if it_end > Parent.npts
                fprintf('Warning***: it_end > length(P)\n')
                it_end = Parent.npts;
            end
            
            obj.Parent=Waveform(Parent.u(1:it_end),Parent.t(1),Parent.dt);
            obj.Daughter=Waveform(Daughter.u(1:it_end),Daughter.t(1),Daughter.dt);
            
            %disp(max(P(1:it_end)));
            
            [obj.time, obj.rf] = ETMTM(...
                obj.Parent.u',obj.Daughter.u',...
                TB,NT,'data',dt,win_len,Poverlap);       
        end
        
        function plot(obj)
            clf;
            subplot(2,1,1)
            plot(obj.Parent); hold on
            plot(obj.Daughter);
            subplot(2,1,2)
            plot(obj.time,obj.rf); hold on
        end
        
    end
    
    methods (Static)
        function test
            fun=zeros(2000,1);
            fun(1350)=1;
            fun=conv(fun,triang(35),'same');
            fun=conv(fun,triang(35),'same');
            fun=conv(fun,triang(35),'same');
            Parent=Waveform(fun,0,0.2);
            fun=zeros(2000,1);
            for ispike = [800, 900, 1000, 1100, 1200, 1300]
                fun(ispike)=randn;
            end
            fun=conv(fun,triang(35),'same');
            fun=conv(fun,triang(35),'same');
            fun=conv(fun,triang(35),'same');
            Daughter=Waveform(fun,0,0.2);
                       
            plot(ReceiverFunction(Parent,Daughter))
                     
        end
        
    end
    
    
end