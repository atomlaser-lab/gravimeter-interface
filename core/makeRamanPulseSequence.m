function makeRamanPulseSequence(dds, varargin)
    %% Default parameters
    delta = 0; % Detuning default value
    t0 = 1e-3; % Default time first pulse
    width = 10e-6; % Default pulse duration
    appliedPhase = 0; % Default applied phase
    dt = 1e-6; % Default time step

    %% Parse input arguments
    if mod(numel(varargin),2) ~= 0
        error('Arguments must appear as name/value pairs!');
    else
        for nn = 1:2:numel(varargin)
            v = varargin{nn+1};
            switch lower(varargin{nn})
                case 'delta'
                    delta = v;
                case 't0'
                    t0 = v;
                case 'width'
                    width = v;
                case {'appliedphase', 'phase'}
                    appliedPhase = v;
                case 'dt'
                    dt = v;
                otherwise
                    error('Option %s not supported', varargin{nn});
            end
        end
    end

    %% Create pulse sequence
    tPulse = (-5*width:dt:5*width)';
    t = tPulse + t0;
    t = t(:); % Ensure t is a column vector

       %% Set powers, phases, and frequencies
    [P, ph, freq] = deal(zeros(numel(t),2));
    
    % Set powers
    pulseStart = t0 - width/2;
    pulseEnd = t0 + width/2;
    P(t >= pulseStart & t <= pulseEnd, :) = 1; % Set to 1 (or the desired power level) during the pulse

    % Set phases
    ph(:,1) = appliedPhase;
    ph(:,2) = appliedPhase;

    % Set frequencies
    freq(:,1) = DDSChannel.DEFAULT_FREQ + delta/4;
    freq(:,2) = DDSChannel.DEFAULT_FREQ - delta/4;


    %% Populate DDS values
    for nn = 1:numel(dds)
        dds(nn).after(t, freq(:,nn), P(:,nn), ph(:,nn));
    end
end
