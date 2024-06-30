function MakeCompositePulse_Rhys(dds,varargin)

%% Set up default variables and parse inputs
t0 = 10e-6;
appliedPhase = 0;
dt = 1e-6;

delta = 0;
PulseType = 'Primitive';
w0 = 2.5e-3;
NoiseAmp = 0;
RampAmp = 0;
NoiseType = 'acceleration';

% % % % % % % New inputs
P1_max = 4;
P2_max = 25;
P_pi = 7;
P_rat = 7;

if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            % % % % Timing
            case 't0'
                t0 = v;
            case 'dt'
                dt = v;
            case 'width'
                tau_seg = v;
            % % % % Beam powers
            case 'p1_max'
                P1_max = v;
            case 'p2_max'
                P2_max = v;
            case 'p_pi'
                P_pi = v;
            case 'p_rat'
                P_rat = v;                

            % % % % Ramps + Noise
            case 'noiseamp'
                NoiseAmp = v;
            case 'rampamp'
                RampAmp = v;
            case 'noisetype'
                NoiseType = v;

            % % % % Light Properties
            case 'w0'
                w0 = v;              

            % % % % Pulse Parameters
            case {'appliedphase','phase'}
                appliedPhase = v;
            case 'delta'
                delta = v;             

            case 'pulsetype'
                PulseType = v;


            % % % % Error Check    
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end



% % % % % Intermediate value
if P2_max*(1 + 1/7) < P1_max*8
    P_tot_max = P2_max*(1 + 1/7);
else
    P_tot_max = P1_max*8;
end


% % % % Pulse Types
if strcmpi(PulseType,'Primitive')
    PulseAreas = [180];
    Phase = [0];
elseif strcmpi(PulseType,'Corpse')
    PulseAreas = [60 300 420];
    Phase = [0 180 0];
elseif strcmpi(PulseType,'BB1')
    PulseAreas = [180 360 180 180];
    Phase = [104.5 313.4 104.5 0];
elseif strcmpi(PulseType,'Knill')
    PulseAreas = [180 180 180 180 180];
    Phase = [240 210 300 210 240];
elseif strcmpi(PulseType,'Waltz')
    PulseAreas = [90 180 270];
    Phase = [0 180 0];
elseif strcmpi(PulseType,'N_90_360_90')
    PulseAreas = [90 360 90];
    Phase = [0 120 0];
elseif strcmpi(PulseType,'Scrofulous')
    PulseAreas = [180 180 180];
    Phase = [60 300 60];
elseif strcmpi(PulseType,'levitt')
    PulseAreas = [90 180 90];
    Phase = [90 0 90];
elseif strcmpi(PulseType,'N_90_240_90')
    PulseAreas = [90 240 90];
    Phase = [240 330 240];
else
    error('Not a listed pulse type')
end
% % % % 


% Two ways to have a desired pulse area: increase power and hold pulse 
% duration constant or increase pulse duration and hold power constant 
% I'll do the former 

% % PulseDurations = round(PulseAreas/180*tau_pi*1e6,0);
numPulseSegments = numel(PulseAreas);

% Convert Pulse Area into pulse power needed
P_tot_segment = PulseAreas/180*P_pi;
P1_segment = P_tot_segment/(1+P_rat);
P2_segment = P_tot_segment/(1+1/P_rat);

if any(P1_segment > P1_max)
    error('Beam 1 doesnt have enough power for a segment')
end
if any(P2_segment > P2_max)
    error('Beam 2 doesnt have enough power for a segment')
end

% Convert beam power into AOM setting
P_AOM1 = P1_segment/P1_max;
P_AOM2 = P2_segment/P2_max;


%% Create Time Vector
% % % Check for DDS Errors
if dt < 1e-6
    error('DDS error: instructions duration less than 1 us')
elseif mod(tau_seg,1e-6) < 1e-6 && mod(tau_seg,1e-6) ~= 0
    error('DDS error: Pi-Pulse duration requires DDS instruction less than 1 us')    
end

OffInstructionDuration = 1e-6;
PulseDuration = numPulseSegments*tau_seg;
tPulse = ((-4*dt): dt : (PulseDuration + 4*dt))';
tPulse = [tPulse + OffInstructionDuration; tPulse(end) + 2*OffInstructionDuration];
t = tPulse + t0;
%% Create Pulses
%
% Set powers, phases, and frequencies
%
[P,ph,freq] = deal(zeros(numel(t),2));
for nn = 1:numPulseSegments
    if numPulseSegments == 1
        t_RisingEdge = t0 + (nn-1)*tau_seg ;%+ (nn ~= 1)*dt
        t_FallingEdge = t0 + nn*tau_seg;  
    else
        t_RisingEdge = t0 + (nn-1)*tau_seg ;%+ (nn ~= 1)*dt
        t_FallingEdge = t0 + nn*tau_seg - dt;        
    end
    %
    % Set powers
    %
    PulseStart = find(abs(t - (t_RisingEdge)) < dt/2,1,'last');
    PulseEnd = find(abs(t - (t_FallingEdge)) < dt/2,1,'first');

    SquareShape = zeros(length(t),1);
    SquareShape(PulseStart:PulseEnd) = 1;
    idx = SquareShape == 1;

    P(:,2) = P(:,2) + P_AOM2(nn).*SquareShape;
    P(:,1) = P(:,1) + P_AOM1(nn).*SquareShape;

%     plot(t,P);
    %
    % Set phases
    %
    ph(idx,2) = Phase(nn);

    %
    % Set frequencies. Channel 1 corresponds to the sideband, channel 2 to
    % the carrier. AOMs are aligned for the +1 order.  If using the +1
    % sideband, and if the R&S synthesizer is set to HFS - DF, then DELTA =
    % DF*2;
    %
    freq(:,1) = DDSChannel.DEFAULT_FREQ + delta/4;
    freq(:,2) = DDSChannel.DEFAULT_FREQ - delta/4;   
end

freq(freq == 0) = DDSChannel.DEFAULT_FREQ;

%% Populate DDS values
for nn = 1:numel(dds)
    if nn == 1
        dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));        
    elseif nn == 2
        dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
    end
end

end