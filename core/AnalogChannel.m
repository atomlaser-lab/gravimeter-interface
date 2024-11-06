classdef AnalogChannel < TimingControllerChannel
    %AnalogChannel Defines an analog channel as an extension
    %to the TimingControllerChannel class.
    %
    %At the moment, there is nothing to extend...

    methods
        function ch = AnalogChannel
            ch = ch@TimingControllerChannel;
            ch.setBounds([-Inf,Inf]);
            ch.IS_ANALOG = true;
        end
    end




end