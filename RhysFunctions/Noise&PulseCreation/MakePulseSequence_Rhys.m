function MakePulseSequence_Rhys(dds,varargin)

%% Set up default variables and parse inputs
f = 384.224e12;
k = 2*pi*f/const.c;

t0 = 10e-6;
width = 30e-6;
width2 = 30e-6;
T = 1e-3;
appliedPhase = 0;
power1 = 0.05*[1,2,1];
dt = 1e-6;
chirp = 2*k*9.795/(2*pi);

delta = 0;
PulseShape = 'square';
w0 = 2.5e-3;
NoiseAmp = 0;
RampAmp = 0;
NoiseType = 'acceleration';

if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            % % % % Timing
            case 't0'
                t0 = v;
            case 't'
                T = v;
            case 'dt'
                dt = v;
            case 'width'
                width = v;
            case 'width2'
                width2 = v;
            % % % % Ramps + Noise
            case 'noiseamp'
                NoiseAmp = v;
            case 'rampamp'
                RampAmp = v;
            case 'noisetype'
                NoiseType = v;
            case 'chirp'
                chirp = v;

            % % % % Light Properties
            case 'w0'
                w0 = v;               
            case 'power1'
                power1 = v;
                if any(power1 < 0)
                    error('Power needs a value between 0 and 1.');
                elseif power1 > 1
                    error('Power needs a value between 0 and 1.');
                else
                    power1 = v;
                end
            case 'power2'
                power2 = v;
                if any(power2 < 0)
                    error('Power needs a value between 0 and 1.');
                elseif power2 > 1
                    error('Power needs a value between 0 and 1.');
                else
                    power2 = v;
                end   

            % % % % Pulse Parameters
            case {'appliedphase','phase'}
                appliedPhase = v;
            case 'delta'
                delta = v;             

            case 'pulseshape'
                if strcmpi(v,'gaussian') == 0 && strcmpi(v,'square') == 0
                    error('Pulse Type must be "Square" or "Gaussian"')
                end
                PulseShape = v;


            % % % % Error Check    
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end


%% Calculate intermediate values
numPulses = max(sum(power2 ~= 0),sum(power1 ~= 0)); % need max incase ch1 or ch2 is off

if numel(appliedPhase) == 0
    tmp = zeros(1,numPulses);
    tmp(end) = appliedPhase;
    appliedPhase = tmp;
end

%% Create Time Vector
% % % Check for DDS Errors
if dt < 1e-6
    error('DDS error: instructions duration less than 1 us')
elseif round(mod(T,1e-6),7) < 1e-6 && round(mod(T,1e-6),7)~= 0
    error('DDS error: Pulse separation time requires DDS instruction less than 1 us')
elseif mod(width,1e-6) < 1e-6 && mod(width,1e-6) ~= 0
    error('DDS error: Pulse duration requires DDS instruction less than 1 us')    
end

OffInstructionDuration = 1e-6;
tPulse = ((-4*dt): dt : (width + 4*dt))';
tPulse = [tPulse + OffInstructionDuration; tPulse(end) + 2*OffInstructionDuration];

tPulse2 = ((-4*dt): dt : (width2 + 4*dt))';
tPulse2 = [tPulse2 + OffInstructionDuration; tPulse2(end) + 2*OffInstructionDuration];

% t = repmat(tPulse,1,numPulses);
t1 = tPulse + t0;
t2 = tPulse2 + t0 + width + T;
t3 = tPulse + t0 + width + T + width2 + T;

t = [t1; t2; t3];

% for  nn = 1:numPulses
%     if numPulses == 1
%         % If ch1 or ch2 is off while the other is on, use the on channel to
%         % determine which pulse should be on
%         if sum(power1 ~= 0) ~= 0
%             nn = find(power1 ~= 0);
%         elseif sum(power1 ~= 0) == 0
%             nn = find(power2 ~= 0);
%         else
%             warning('no power during ya pulses')
%         end
%         t(:,1) = t(:,1) + t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width;
%     else
%         t(:,nn) = t(:,nn) + t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width;
%     end
% end
% t = t(:);

%% Create Noise 
% % % % Create intensitity noise profile
if strcmpi(NoiseType,'acceleration')
    r = 0.5*NoiseAmp*t.^2;
    I_Noise_factor = exp(-2*r.^2/w0^2);
elseif strcmpi(NoiseType,'white')
    I_Noise_factor = NoiseAmp*normrnd(0,1,length(t),1);
elseif strcmpi(NoiseType,'white') == 0 && strcmpi(NoiseType,'acceleration')
    warning('Type must be "white" or "acceleration"')
    return
end

%% Create Ramp
% % % % Create intensitity ramp
r = 0.5*RampAmp*t.^2;
I_Ramp_factor = exp(r.^2/w0^2);

%% Combine Noise and Ramps
I_factor = I_Ramp_factor.*I_Noise_factor;

%% Create Pulses
%
% Set powers, phases, and frequencies
%
[P,ph,freq] = deal(zeros(numel(t),2));
for nn = 1:numPulses
    if numPulses == 1
        if sum(power1 ~= 0) ~= 0
            nn = find(power1 ~= 0);
        elseif sum(power1 ~= 0) == 0
            nn = find(power2 ~= 0);
        else
            warning('no power during ya pulses')
        end
        t_RisingEdge = t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width2;
        t_FallingEdge = t0 + width + (nn-1)*T + (nn~=1)*width2 + (nn==3)*width;
    else
        t_RisingEdge = t0 + (nn-1)*T + (nn~=1)*width + (nn==3)*width2;
        t_FallingEdge = t0 + width + (nn-1)*T + (nn~=1)*width2 + (nn==3)*width;
    end
    %
    % Set powers
    %
    if strcmpi(PulseShape,'Gaussian') == 1
        fwhm = width/(2*sqrt(log(2)));
        P(:,1) = P(:,1) + power1(nn).*I_factor(:,1).*exp(-(t - t_RisingEdge).^2/fwhm.^2);
        P(:,2) = P(:,2) + power2(nn).*I_factor(:,1).*exp(-(t - t_RisingEdge).^2/fwhm.^2);
    elseif strcmpi(PulseShape,'Square') == 1
        PulseStart = find(abs(t - (t_RisingEdge + OffInstructionDuration)) < dt/2,1,'last');
        PulseEnd = find(abs(t - (t_FallingEdge + OffInstructionDuration)) < dt/2,1,'first');        
        SquareShape = zeros(length(t),1);
        SquareShape(PulseStart:PulseEnd) = 1;
        idx = SquareShape == 1;
        P(:,2) = P(:,2) + power2(nn).*I_factor(:,1).*SquareShape;
        P(:,1) = P(:,1) + power1(nn).*I_factor(:,1).*SquareShape;
    end
    %
    % Set phases
    %
   ph(idx,2) = appliedPhase(nn);


    %
    % Set frequencies
    %
    freq(:,1) = DDSChannel.DEFAULT_FREQ + delta/4 + 0.25/1e6*(chirp*t);
    freq(:,2) = DDSChannel.DEFAULT_FREQ - delta/4 - 0.25/1e6*(chirp*t);   
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