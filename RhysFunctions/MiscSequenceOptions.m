classdef MiscSequenceOptions < SequenceOptionsAbstract
    
    properties
        tof2
        detuning2
        P1Max
        P2Max
        AbsAnalysis_AllROI
        DropCamera
    end

    methods
        function self = MiscSequenceOptions(varargin)
            self.setDefaults;
            self = self.set(varargin{:});
        end

        function self = setDefaults(self)
            self.tof2 = 7e-3;
            self.detuning2 = 0;
            self.P1Max = 5;
            self.P2Max = 5;
            self.AbsAnalysis_AllROI = 1;
            self.DropCamera = 'in-trap';
        end
    end
end