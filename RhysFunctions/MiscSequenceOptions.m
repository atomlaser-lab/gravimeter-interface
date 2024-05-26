classdef MiscSequenceOptions < SequenceOptionsAbstract
    
    properties
        tof2
        detuning2
    end

    methods
        function self = MiscSequenceOptions(varargin)
            self.setDefaults;
            self = self.set(varargin{:});
        end

        function self = setDefaults(self)
            self.tof2 = 7e-3;
            self.detuning2 = 0;
        end
    end
end