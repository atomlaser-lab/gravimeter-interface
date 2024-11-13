classdef AnalogChannel < TimingControllerChannel
    %AnalogChannel Defines an analog channel as an extension
    %to the TimingControllerChannel class.
    %
    %At the moment, there is nothing to extend...

    properties(Constant)
        VOLTAGE_BOUNDS = [-10,10];
    end

    methods
        function ch = AnalogChannel
            ch = ch@TimingControllerChannel;
            ch.setBounds([-Inf,Inf]);
            ch.IS_ANALOG = true;
        end

        function values_out = convert(ch,values_in)
            values_out = convert@TimingControllerChannel(ch,values_in);
            if any(values_out < AnalogChannel.VOLTAGE_BOUNDS(1) | values_out > AnalogChannel.VOLTAGE_BOUNDS(2))
                error('Converted voltage values are outside the range [%.0f,%.0f] V!',AnalogChannel.VOLTAGE_BOUNDS(1),AnalogChannel.VOLTAGE_BOUNDS(2));
            end
        end
    end




end