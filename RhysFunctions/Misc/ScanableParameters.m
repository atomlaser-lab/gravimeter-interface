classdef ScanableParameters < handle
    properties(Constant)
        Voltage = 'Voltage (V)';
        TOF = 'TOF (ms)';
        PulseDuration = '\tau (us)';
        TwoPhoton = '\delta (kHz)';        
        Power = 'Total Power (mW)';
        Run = 'Run Number (Arb)';
        Phase = 'DDS Phase (deg)'; 
        T = 'T (ms)';
    end

    methods(Static)

    end
end