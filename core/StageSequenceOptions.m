classdef StageSequenceOptions < SequenceOptionsAbstract
    properties
        mot
        cmot
        pgc
        pump
        mag
        evap_mag
        dipoles
        evap_dipoles

        override
    end

    methods
        function self = StageSequenceOptions(varargin)
            self.setDefaults;
            self.set(varargin{:});
        end

        function self = setDefaults(self)
            self.override = 0;
            self.enable_all_stages;
        end

        function self = enable_all_stages(self)
            self.mot = 1;
            self.cmot = 1;
            self.pgc = 1;
            self.pump = 1;
            self.mag = 1;
            self.evap_mag = 1;
            self.dipoles = 1;
            self.evap_dipoles = 1;
        end

        function self = disable_all_stages(self)
            self.mot = 0;
            self.cmot = 0;
            self.pgc = 0;
            self.pump = 0;
            self.mag = 0;
            self.evap_mag = 0;
            self.dipoles = 0;
            self.evap_dipoles = 0;
        end

        function self = disable_dipoles(self)
            self.dipoles = 0;
            self.evap_dipoles = 0;
        end

        function self = enable_only_laser_cooling(self)
            self.disable_all_stages;
            self.mot = 1;
            self.cmot = 1;
            self.pgc = 1;
        end

        function r = use_mot(self)
            r = self.mot;
        end

        function r = use_cmot(self)
            r = self.cmot & (self.use_mot | self.override);
        end

        function r = use_pgc(self)
            r = self.pgc & (self.use_cmot | self.override);
        end

        function r = use_pump(self)
            r = self.pump & (self.use_pgc | self.override);
        end

        function r = use_mag(self)
            r = self.mag & (self.use_pgc | self.override);
        end

        function r = use_evap_mag(self)
            r = self.evap_mag & (self.use_mag | self.override);
        end

        function r = use_dipoles(self)
            r = self.dipoles & (self.evap_mag | self.override);
        end

        function r = use_evap_dipoles(self)
            r = self.evap_dipoles & (self.use_dipoles | self.override);
        end
    end
end