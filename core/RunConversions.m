classdef RunConversions < handle
    methods(Static)


        function V = UD_A_to_V(I)
            if I > 16
                warning('UD bias current limit is 17 A. Current has been set to 17 A')
                I = 17;
            end

            if I <= 13.2
                V = 0.166667 * I + 0.000000;
            else
                V = 0.163366 * I + 0.041584;
            end
        end


        function V = mot_power(power)
            % This is a percentage of total power
            Voltage = [0,1.4000,1.5000,1.6000,1.7000,1.8000,1.9000,2.0000,2.2500,2.5000,3.0000,3.5000,4.0000,4.5000,4.7500,5.0000];
            NormalisedPower = [0,0.000002889447236,0.000056532663317,0.000383165829146,0.001482412060302,0.004020100502513,0.008793969849246,0.015954773869347,0.051256281407035,0.109296482412060,0.300125628140704,0.542713567839196,0.747487437185930,0.917085427135678,0.967336683417085,1.000000000000000];
            func = @(DesiredPower) interp1(NormalisedPower,Voltage,DesiredPower);
            if power > 1 || power < 0
                error('MOT power must be between zero and one')
            end
            V = func(power);
            if power == 0
                V = 0;
            end
        end

        function V = repump_power(power)
            % This is a percentage of total power
            Voltage = [1.2000,1.3000,1.5000,2.0000,2.5000,3.0000,3.5000,4.0000,4.5000,5.0000,5.5000,6.0000,6.5000,7.0000,7.5000,8.0000,8.5000];
            NormalisedPower = [0,0.000042253521127,0.002521126760563,0.046478873239437,0.129577464788732,0.225352112676056,0.330985915492958,0.436619718309859,0.535211267605634,0.619718309859155,0.718309859154930,0.788732394366197,0.859154929577465,0.901408450704225,0.957746478873239,0.971830985915493,1.000000000000000];
            func = @(DesiredPower) interp1(NormalisedPower,Voltage,DesiredPower);
            if power > 1 || power < 0
                error('MOT power must be between zero and one')
            end
            V = func(power);
            if power == 0
                V = 0;
            end
        end

        function V = dipole25(power)
            V = (power + 0.15)/2.76;    % 12/04/2024, with servo
            %             V = 0.325898426677744*power + 1.127093693591155; % 18/12/2023
        end

        function V = dipole50(power)
            V = (power + 0.19)/5.42;    % 12/04/2024, with servo
            %             V =  0.167016409500481*power + 0.875247148388929; % 18/12/2023
        end

        function V = mot_coil(current)
            %             V = (current - 0.3)/6;
            %             V(current == 0) = -0.075;
            V = (current - 0.45)/6;
        end

        function I = mot_coil_reverse(V)
            %             V = (current - 0.3)/6;
            %             V(current == 0) = -0.075;
            I = 6*V + 0.45;
        end

        function V = imaging(detuning)
            % % %             V = -detuning*0.472/6.065 + 8.533;
            % %             func = @(x) (82.6202 - 8.7259*x + 2.0478*x.^2 - 0.0827*x.^3);
            % % %             f = (detuning + 211.79)/2;
            % %             f = (detuning + 106.4645*2)/2;
            % % %             func = @(x) (8.7236-.8121*x+0.0075*x.^2 -0.0031*x.^3);
            % % %             f = (detuning + 105.9828)*2;
            % %             xx = linspace(0,10,101);
            % %             V = interp1(func(xx),xx,f,'pchip');
            frequency = [93.248299, 97.059586, 98.692995, 99.5097, 100.326404, 101.143109, ...
                102.232048, 103.048752, 103.865457, 104.682161, 105.311338, ...
                106.162672, 107.014005, 107.865338, 108.716671, 109.357566, ...
                110.295859, 110.966067, 111.770318, 112.574569, 113.378819, ...
                114.18307, 114.987321, 115.791571, 116.595822, 117.400073];

            voltage = [7, 7.5, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, ...
                8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10];

            % Interpolating function
            func = @(f) interp1(frequency, voltage, f, 'pchip');
            Abs_frequency = (detuning-.4+106.4645*2)/2;
            V = func(Abs_frequency);
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