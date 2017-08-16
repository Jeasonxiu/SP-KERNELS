classdef nth_deriv
    properties
        n
        fx
        x
        fx_deriv
    end

    methods
        function obj = nth_deriv(n,fx,x)
            obj.n=n;
            obj.fx=fx;
            obj.x=x;
            obj = compute_deriv(obj);
        end
        function obj = compute_deriv(obj)
            FX=fftshift(fft(obj.fx));
            dx=obj.x(2)-obj.x(1);
            freq_max=1/2/dx;
            f=linspace(-freq_max,freq_max,length(FX));
            FX_DERIV=FX.*(2*pi*1i*f).^obj.n;
            obj.fx_deriv=ifft(ifftshift(FX_DERIV));
            obj.fx_deriv=real(obj.fx_deriv);
        end
        
        function plot(obj)
            figure(1)
            clf;
            plot(obj.x,obj.fx); hold on;
            plot(obj.x,obj.fx_deriv)
            label1='f(x)';
            label2=sprintf('f^{(%.1f)}(x)',obj.n);
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
            fx=exp(-x.^2);
            c=nth_deriv(n,fx,x);
            plot(c);
            fx_check=-2*exp(-x.^2).*x;
            plot(c.x,fx_check,'--k');
        end
    end
end