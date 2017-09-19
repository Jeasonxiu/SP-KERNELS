classdef myhilbert
    properties
        fx
        x
        fx_hi
        env
    end

    methods
        function obj = myhilbert(fx,x)
            obj.fx=fx;
            obj.x=x;
            obj = compute_hilbert(obj);
        end
        function obj = compute_hilbert(obj)
            FX=fftshift(fft(obj.fx));
            dx=obj.x(2)-obj.x(1);
            freq_max=1/2/dx;
            f=linspace(-freq_max,freq_max,length(FX));
            FX_HI=FX*1i.*sign(f);
            obj.fx_hi=ifft(ifftshift(FX_HI));
            obj.fx_hi=real(obj.fx_hi);
            obj.env=sqrt(obj.fx.^2+obj.fx_hi.^2);
        end
        
        function plot(obj)
            figure(1)
            clf;
            plot(obj.x,obj.fx); hold on;
            plot(obj.x,obj.fx_hi)
            plot(obj.x,obj.env,'--k')
            label1='f(x)';
            label2=sprintf('f_{hi}(x)');
            legend(label1,label2)
            xlabel('x')
        end
    end
    
    methods (Static)
        function test()
            npts=2000;
            n=1.0;
            x=linspace(-10,10,npts);
            %fx=x;
            fx=exp(-x.^2).*cos(3*x);
            c=nth_deriv(2,fx,x);
            c=myhilbert(c.fx_deriv,c.x);
            plot(c);
            %fx_check=-2*exp(-x.^2).*x;
            %plot(c.x,fx_check,'--k');
        end
    end
end