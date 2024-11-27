classdef SequenceOptions < SequenceOptionsAbstract
    %SEQUENCEOPTIONS Defines a class for passing options to the
    %gravimeter make sequence function
    
    properties
        %
        % Preparation properties
        %
        load_time
        detuning
        redpower
        raycus
        tof
        %
        stage
        % Other properties
        %
        nd
        dkc
        params
        %%% JM ADDED THIS
        param1
        param2
        param3
    end
    
    methods
        function self = SequenceOptions(varargin)
            self.setDefaults;
            self = self.set(varargin{:});
        end

        function self = setDefaults(self)
            self.load_time = 7.5;
            self.detuning = 0;
            self.redpower = 2;
            self.raycus = 2;
            self.tof = 20e-3;
            self.params = [];
            %%% JM ADDED THIS
            self.param1 = [];
            self.param2 = [];
            self.param3 = [];
           
            self.stage = StageSequenceOptions;
            self.nd = FeedbackOptions;
            self.dkc = DKC_Options;
        end
        
        function self = set(self,varargin)
            set@SequenceOptionsAbstract(self,varargin{:});

            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            else
                for nn = 1:2:numel(varargin)
                    switch lower(varargin{nn})
                        case 'dipoles'
                            self.raycus = varargin{nn+1};
                            self.redpower = varargin{nn+1};
                    end
                end
            end
        end

        
    end
    
end