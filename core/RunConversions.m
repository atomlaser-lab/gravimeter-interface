classdef RunConversions < handle
    methods(Static)
        function V = dipole25(power)
            V = (power + 0.1795)/2.7065;
        end
        
        function V = dipole50(power)
            V = (power + 66.9e-3)/4.9909;
        end
        
        function V = mot_coil(current)
%             V = (current - 0.3)/6;
%             V(current == 0) = -0.075;
            V = (current - 0.45)/6;
        end
        
        function V = imaging(detuning)
%             V = -detuning*0.472/6.065 + 8.533;
            func = @(x) (82.6202 - 8.7259*x + 2.0478*x.^2 - 0.0827*x.^3);
            f = (detuning + 211.79)/2;

%             func = @(x) (8.7236-.8121*x+0.0075*x.^2 -0.0031*x.^3);
%             f = (detuning + 105.9828)*2;
            xx = linspace(0,10,101);
            V = interp1(func(xx),xx,f,'pchip');
        end
        
        function V = microwave(detuning)
            V = 7.9141 - 0.0152*detuning - 2.555e-5*detuning.^2 + 3.6751e-7*detuning.^3;
        end
        
        function V = mot_freq(detuning)
            func = @(x) (53.051 + 8.6164*x - 1.5183*x.^2 + 0.24203*x.^3 - 0.010976*x.^4);
            f = (detuning + 211.79)/2;
            xx = linspace(0,10,101);
            V = interp1(func(xx),xx,f,'pchip');
        end
        
        function V = repump_freq(detuning)
            func = @(x) 51.933 + 9.3739*x - 1.8124*x.^2 + 0.28129*x.^3 - 0.01269*x.^4;
            f = detuning + 78.47;
            xx = linspace(0,10,101);
            V = interp1(func(xx),xx,f,'pchip');
        end
    end
end