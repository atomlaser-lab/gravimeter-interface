classdef DKC_Options < SequenceOptionsAbstract
    properties
        type
        power
        delay
        duration
        dipoles
    end

    methods
        function self = DKC_Options(varargin)
            self.setDefaults;
            self.set(varargin{:});
        end

        function self = setDefaults(self)
            self.type = 'guiding';
            self.power = 11;
            self.delay = 0;
            self.duration = 20e-3;
            self.dipoles = [];
        end

        function r = is_guide(self)
            if any(strcmpi(self.type,{'guide','guiding'}))
                r = 1;
            else
                r = 0;
            end
        end

        function r = is_kick(self)
            if any(strcmpi(self.type,{'kick','normal'}))
                r = 1;
            else
                r = 0;
            end
        end

        function r = is_test(self)
            if any(strcmpi(self.type,{'test','testing'}))
                r = 1;
            else
                r = 0;
            end
        end
    end
end