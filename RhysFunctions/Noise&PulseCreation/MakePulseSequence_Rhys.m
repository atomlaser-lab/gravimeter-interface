function MakePulseSequence_Rhys(dds,varargin)

%% Set up default variables and parse inputs
f = 384.224e12;
k = 2*pi*f/const.c;
t0 = 10e-3;
width = 30e-6;
T = 1e-3;
appliedPhase = 0;
power1 = 0.05*[1,2,1];
chirp = 2*k*9.795/(2*pi);
order = 1;
start_order = 0;
dt = 1e-6;

delta = 0;
PulseType = 'Gaussian';
w0 = 12.5e-3; %("~2 cm FWHM" section 4.2 of Hardman)
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
            case 'f'
                f = v;
                k = 2*pi*f/const.c;
            case 'k'
                k = v;                
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
            case 'chirp'
                chirp = v;

            case 'pulsetype'
                if strcmpi(v,'gaussian') == 0 && strcmpi(v,'square') == 0
                    error('Pulse Type must be "Square" or "Gaussian"')
                end
                PulseType = v;


            % % % % Error Check    
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end


%% Calculate intermediate values
recoil = const.hbar*k^2/(2*const.mRb*2*pi);
Stark = 0.1304;


numPulses = max(sum(power2 ~= 0),sum(power1 ~= 0));
fwhm = width/(2*sqrt(log(2)));

if numel(appliedPhase) == 0
    tmp = zeros(1,numPulses);
    tmp(end) = appliedPhase;
    appliedPhase = tmp;
end

%% Checks
if T < width && numPulses>1
    error('Pulse separation time is less than the pulse duration: pulses are not separated')
end
if t0 < width
%     warning('initial drop time less than pulse width: pulse starts before the desired initial drop time')
end


%% Create Time Vector
% % % Check for DDS Errors
if dt < 1e-6
    error('DDS error: instructions duration less than 1 us')
elseif mod(T,1e-6) < 1e-6 && mod(T,1e-6)~= 0
    error('DDS error: Pulse separation time requires DDS instruction less than 1 us')
elseif mod(width,1e-6) < 1e-6 && mod(width,1e-6) ~= 0
    error('DDS error: Pulse duration requires DDS instruction less than 1 us')    
end

OffInstructionDuration = 1e-6;
tPulse = (0: dt :width)';
tPulse = [0;  tPulse + OffInstructionDuration; tPulse(end) + 2*OffInstructionDuration];
t = repmat(tPulse,1,numPulses);

% for  nn = 1:numPulses
%     t(:,nn) = t(:,nn) + t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width;
% end

for  nn = 1:numPulses
    if numPulses == 1
        nn = find(power1 ~= 0);
        t(:,1) = t(:,1) + t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width;
    else
        t(:,nn) = t(:,nn) + t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width;
    end
end

t = t(:);

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
I_Ramp_factor = exp(2*r.^2/w0^2);

%% Combine Noise and Ramps
I_factor = I_Ramp_factor.*I_Noise_factor;

%% Create Pulses
%
% Set powers, phases, and frequencies
%
[P,ph,freq] = deal(zeros(numel(t),2));
for nn = 1:numPulses
    if numPulses == 1
        nn = find(power1 ~= 0);
        t_RisingEdge = t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width;
    else
        t_RisingEdge = t0 + (nn-1)*T + (nn~=3)*(nn-1)*width + (nn==3)*(nn-1)*width;
    end
    %
    % Set powers
    %
    if strcmpi(PulseType,'Gaussian') == 1
        P(:,1) = P(:,1) + power1(nn).*I_factor(:,1).*exp(-(t - t_RisingEdge).^2/fwhm.^2);
        P(:,2) = P(:,2) + power2(nn).*I_factor(:,1).*exp(-(t - t_RisingEdge).^2/fwhm.^2);
    elseif strcmpi(PulseType,'Square') == 1
        PulseStart = find(abs(t - (t_RisingEdge + OffInstructionDuration)) < dt/2,1,'last');
        PulseEnd = find(abs(t - (t_RisingEdge + width + OffInstructionDuration)) < dt/2,1,'first');        
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
    % Set frequencies. Channel 1 corresponds to the sideband, channel 2 to
    % the carrier. AOMs are aligned for the +1 order.  If using the +1
    % sideband, and if the R&S synthesizer is set to HFS - DF, then DELTA =
    % DF*2;
    %
    freq(idx,1) = DDSChannel.DEFAULT_FREQ + delta/4;
    freq(idx,2) = DDSChannel.DEFAULT_FREQ - delta/4;   
end

freq(freq == 0) = DDSChannel.DEFAULT_FREQ;
% freq(freq ~= 0) = 50;

%% Populate DDS values
for nn = 1:numel(dds)
    if nn == 1
        dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));        
    elseif nn == 2
%         dds(nn).after(t+1e-6,freq(:,nn),P(:,nn),ph(:,nn));
        dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
    end
end


% Require two off values to actually turn the pulse off
% dds(1).after(1e-6,110,0,0);
% dds(2).after(1e-6,110,0,0);

end