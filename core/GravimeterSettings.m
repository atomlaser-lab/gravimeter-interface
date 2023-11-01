classdef GravimeterSettings
    %GRAVIMETEROPTIONS Defines a class for passing options to the
    %gravimeter make sequence function
    
    properties
        %
        % Preparation properties
        %   
        
        dipole_final_value
        motCoilOff
        
       
        %
        % Velocity selection properties
        %
        
        width_vs
        vs_power
        tvs
        
        % Interferometer properties
        %
        T_exp
        T_int
        chirp
        final_phase
        bragg_power
        width_Bragg
        Tasym
        Tsep
        BraggPower_custom
        phase_custom
        
        %
        %Imaging settings
        %
        detuning
        Drop_time
        
        %         OD_max
        %         filter
        
        
        %
        % Other properties
        %
        params
        rc
    end
    
    methods
        function self = GravimeterSettings(varargin)
            self = self.set(varargin{:});
        end
        
        function self = set(self,varargin)
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            else
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        
                        case {'Imaging detuning','detuning'}
                            self.detuning = v;
                        case {'dipole_param','dipole','final_dipole_power'}
                            self.dipole_final_value = v;
                        case 'motcoiloff'
                            self.motCoilOff = v;
                        case {'tof','droptime','drop_time'}
                            self.Drop_time = v;
                        case 't0'
                            self.T_exp = v;
                        case {'t_int','t'}
                            self.T_int = v;
                        case {'final_phase','phase'}
                            self.final_phase = v;
                        case {'bragg_power','power'}
                            self.bragg_power = v;
                        case 'width_bragg' 
                            self.width_Bragg = v;
                        case {'tasym','asym'}
                            self.Tasym = v;
                        case {'tsep','separation'}
                            self.Tsep = v;
                        case 'chirp'
                            self.chirp = v;
                        case 'width_vs'
                            self.width_vs = v;
                        case 'vs_power'
                            self.vs_power = v;
                        case 'tvs'
                            self.tvs = v;
                        case 'braggpower_custom'
                            self.BraggPower_custom = v;
                        case 'phase_custom'
                            self.phase_custom = v;
                        case 'params'
                            self.params = v;
                        case 'rc'
                            self.rc = v;
                        otherwise
                            warning('Option ''%s'' not supported',varargin{nn})
                    end
                end
            end
        end
        
        function self = replace(self,opt)
            p = properties(opt);
            for nn = 1:numel(p)
                if ~isempty(opt.(p{nn}))
                    self.(p{nn}) = opt.(p{nn});
                end
            end
        end
        
        function s = print(self)
            p = properties(self);
            sargs = {};
            for nn = 1:numel(p)
                v = self.(p{nn});
                if isempty(v)
                    sargs{nn} = sprintf('''%s'',%s',p{nn},'[]');
                elseif ischar(v) || isstring(v)
                    sargs{nn} = sprintf('''%s'',''%s''',p{nn},v);
                else
                    sargs{nn} = sprintf('''%s'',%.6g',p{nn},v);
                end
            end
            s = strjoin(sargs,',');
            s = sprintf('GravimeterSettings(%s);',s);
        end
        
    end
    
end