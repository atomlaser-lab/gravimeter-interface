classdef RamanSequenceOptions < SequenceOptionsAbstract
    
    properties
        width
        power
        df
        OnOff
    end

    methods
        function self = RamanSequenceOptions(varargin)
            self.setDefaults;
            self = self.set(varargin{:});
        end

        function self = setDefaults(self)
            self.width = 10e-6;
            self.power = 0.5;
            self.df = 0;
            self.OnOff = 0;
        end
    end
end