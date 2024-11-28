classdef StageSequenceOptions < SequenceOptionsAbstract
    %STAGESEQUENCEOPTIONS Defines a set of properties and methods for
    %turning stages on and off in a sequence maker file
    properties
        mot             %Enable MOT stage
        cmot            %Enable CMOT stage
        pgc             %Enable PGC stage
        pump            %Enable optical pumping stage
        mag             %Enable loading into magnetic trap
        evap_mag        %Enable evaporation in magnetic trap
        dipoles         %Enable dipole trap
        evap_dipoles    %Enable evaporation in optical dipole trap

        override        %Override hierarchy of stages in use_* functions    
    end

    methods
        function self = StageSequenceOptions(varargin)
            %STAGESEQUENCEOPTIONS Constructs a StageSequenceOptions object
            %
            %   SELF = StageSequenceOptions(varargin) constructs
            %   StageSequenceOptions using variable argument list
            %   conforming to the SequenceOptionsAbstract method
            self.setDefaults;
            self.set(varargin{:});
        end

        function self = setDefaults(self)
            %SETDEFAULTS Sets default values
            self.override = 0;
            self.enable_all_stages;
        end

        function self = enable_all_stages(self)
            %ENABLE_ALL_STAGES Enables all stages
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
            %DISABLE_ALL_STAGES Disables all stages
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
            %DISABLE_DIPOLES Disables the dipole and dipole evaporation
            %stages
            self.dipoles = 0;
            self.evap_dipoles = 0;
        end

        function self = enable_only_laser_cooling(self)
            %ENABLE_ONLY_LASER_COOLING Enables only the laser cooling
            %stages, excluding the optical pumping stage
            self.disable_all_stages;
            self.mot = 1;
            self.cmot = 1;
            self.pgc = 1;
        end

        function r = use_mot(self)
            %USE_MOT Checks if MOT is in use
            %
            %   R = USE_MOT() Returns true if MOT is in use, false
            %   otherwise
            r = self.mot;
        end

        function r = use_cmot(self)
            %USE_CMOT Checks if CMOT is in use
            %
            %   R = USE_CMOT() Returns true if MOT and CMOT are in use,
            %   false otherwise
            r = self.cmot & (self.use_mot | self.override);
        end

        function r = use_pgc(self)
            %USE_PGC Checks if PGC is in use
            %
            %   R = USE_PGC() Returns true if MOT, CMOT, and PGC are in
            %   use, false otherwise
            r = self.pgc & (self.use_cmot | self.override);
        end

        function r = use_pump(self)
            %USE_PUMP Checks if optical pumping is in use
            %
            %   R = USE_PUMP() Returns true if all stages up to and
            %   including optical pumping are in use, false otherwise
            r = self.pump & (self.use_pgc | self.override);
        end

        function r = use_mag(self)
            %USE_MAG Checks if the magnetic trap is in use
            %
            %   R = USE_MAG() Returns true if all stages up to and
            %   including the magnetic trap are in use, false otherwise
            r = self.mag & (self.use_pump | self.override);
        end

        function r = use_evap_mag(self)
            %USE_EVAP_MAG Checks if the magnetic trap evaporation is in use
            %
            %   R = USE_MAG() Returns true if all stages up to and
            %   including the magnetic trap evaporation are in use, false
            %   otherwise
            r = self.evap_mag & (self.use_mag | self.override);
        end

        function r = use_dipoles(self)
            %USE_DIPOLES Checks if the optical dipoles are in use
            %
            %   R = USE_DIPOLES() Returns true if all stages up to and
            %   including the dipole stage are in use, false otherwise
            r = self.dipoles & (self.use_evap_mag | self.override);
        end

        function r = use_evap_dipoles(self)
            %USE_EVAP_DIPOLES Checks if the optical dipole evaporation is
            %in use
            %
            %   R = USE_EVAP_DIPOLES() Returns true if all stages up to and
            %   including the dipole evaporation are in use, false
            %   otherwise
            r = self.evap_dipoles & (self.use_dipoles | self.override);
        end
    end
end