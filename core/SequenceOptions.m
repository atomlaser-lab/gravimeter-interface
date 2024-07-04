classdef SequenceOptions < SequenceOptionsAbstract
    %SEQUENCEOPTIONS Defines a class for passing options to the
    %make sequence function

    properties
        %
        % Preparation properties. These are native properties to this set
        % of sequence options
        %
        detuning        %Detuning of imaging light in MHz
        dipoles         %Final power for the two dipole beams in W
        tof             %Time-of-flight in ms
%         imaging_type    %Imaging system to use (drop 1, 2, 3, or 4)
        params          %Additional parameters for optimisation
        extraparams
        MOT_LoadTime    %MOT load time in seconds
        %
        repetition


        MOT_status
        CMOT_status
        PGC_status
        LoadMagTrap_status
        MagEvaporation_status
        LoadOpticalTrap_status
        OpticalEvaporation_status
        BECCompression_status
        MagneticInsensitive_status

        JustMOT
        TwoStateImaging
        % These are sub-groupings of options
        %
        StatePrep
        raman
        bragg
        mw
        misc
    end


    methods
        function self = SequenceOptions(varargin)
            %SEQUENCEOPTIONS Create a SequenceOptions object
            self.bragg = BraggSequenceOptions;
            self.mw = MicrowaveSequenceOptions;
            self.misc = MiscSequenceOptions;

            self.setDefaults;
            self = self.set(varargin{:});
        end

        function self = setDefaults(self)
            %SETDEFAULTS Set default property values
            self.detuning = 0;
            self.dipoles = 1.32;
            self.tof = 25e-3;
%             self.imaging_type = 'drop 2';
            self.params = [];
            self.MOT_LoadTime = 4;
            self.extraparams = [];
            self.repetition = [];
            self.MOT_status =1 ;
            self.CMOT_status =1 ;
            self.PGC_status =1;
            self.LoadMagTrap_status =1 ;
            self.MagEvaporation_status = 1;
            self.LoadOpticalTrap_status = 1;
            self.OpticalEvaporation_status = 1;
            self.BECCompression_status = 0;
            self.MagneticInsensitive_status = 0;
            
            self.JustMOT = 0;
            self.TwoStateImaging = 1;
            self.raman = 0;
            self.StatePrep = 0;
            
            self.bragg.setDefaults;
            self.mw.setDefaults;
            self.misc.setDefaults;
        end

        function self = set(self,varargin)
            %SET Sets the options
            %
            %   SELF = SELF.SET(VARARGIN) Sets properties according to
            %   name/value pairs.  For nested options, use name/value pairs
            %   as 'bragg',{'name',value,'name2',value2,...}
            set@SequenceOptionsAbstract(self,varargin{:});

            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            else
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case 'camera'
                            self.imaging_type = v;
                        case 't0'
                            self.bragg.t0 = v;
                        case 'ti'
                            self.bragg.ti = v;
                        case {'tint','t'}
                            self.bragg.T = v;
                        case 'phase'
                            self.bragg.phase = v;
                        case 'power'
                            self.bragg.power = v;
                        case {'tasym','asym'}
                            self.bragg.Tasym = v;
                        case {'tsep','separation'}
                            self.bragg.Tsep = v;
                        case 'chirp'
                            self.bragg.chirp = v;
                    end
                end
            end
        end

    end

end