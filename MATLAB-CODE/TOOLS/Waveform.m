classdef Waveform
    properties
        t=[]
        u=[]
        npts
        dt
    end
    methods
        function obj=Waveform(u,tstart,dt)
            obj.u=u;
            obj.npts=length(obj.u);
            obj.dt=dt;
            obj.t=linspace(0,obj.npts-1,obj.npts)*obj.dt+tstart; 
            if length(obj.t) ~= length(obj.u)
                error('u and t must be the same length')
            end     
        end
        
        function obj=decimate(obj,n)
            unew=decimate(obj.u,n,'FIR');
            tnew=decimate(obj.t,n,'FIR');
            obj=Waveform(unew,tnew(1),obj.dt*n);
        end
        
        function plot(obj)
            plot(obj.t,obj.u)
        end
    end
    methods (Static)
        function test()
            fun=sind((-10000:10000)/8);
            fun=fun-mean(fun);
            Wf=Waveform(fun,0,10);
            hold on
            plot(Wf);
            plot(decimate(Wf,200))
            
        end
    end
    
end