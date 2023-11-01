function MakePulseSequence_Rhys_new(dds,varargin)

%% Set up variables and parse inputs
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
mirrorSwitch = 1;

NumPulseWidths = 5;
PulseType = 'Gaussian';

% % % new
I_factor = 1;
I_ratio = 1;
% This is I_sideband/I_carrier
% I carrier is DDS Amp 1
% % %

if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 't0'
                t0 = v;
            case 't'
                T = v;
            case 'dt'
                dt = v;
            case 'width'
                width = v;
            case {'appliedphase','phase'}
                appliedPhase = v;
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
            case 'chirp'
                chirp = v;
            case 'f'
                f = v;
                k = 2*pi*f/const.c;
            case 'k'
                k = v;
            case 'numpulsewidths'
                NumPulseWidths = v;
            case 'pulsetype'
                if strcmpi(v,'gaussian') == 0 && strcmpi(v,'square') == 0
                    error('Pulse Type must be "Square" or "Gaussian"')
                end
                PulseType = v;
            case 'i_factor'
                I_factor = v;
            case 'i_ratio'
                I_ratio = v;
            case 'delta'
                delta = v;
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end


%% Calculate intermediate values
recoil = const.hbar*k^2/(2*const.mRb*2*pi);
Stark = 4;

numPulses = sum(power1 ~= 0);
fwhm = width/(2*sqrt(log(2)));

if numel(appliedPhase) == 0
    tmp = zeros(1,numPulses);
    tmp(end) = appliedPhase;
    appliedPhase = tmp;
end

%% Checks
if T < width && numPulses>1
    error('Pulse separation time is less than the pulse duration: pulses occur at the same time')
end
if t0 < width
%     error('nitial drop time less than pulse width: pulse starts before the desired initial drop time')
end


%% Create vectors
tPulse = (-NumPulseWidths*width/2 - dt : dt :NumPulseWidths*width/2 + dt)';
tPulse = (-NumPulseWidths*width/2 : dt :NumPulseWidths*width/2)';
tPulse = [tPulse(1) - 1e-6;  tPulse; tPulse(end) + 1e-6];

% I've add extra points between pulses to explicitly set power to zero
t = repmat(tPulse,1,numPulses);
for  nn = 1:numPulses
    t(:,nn) = t(:,nn) + t0 + (nn-1)*T;
end
t = t(:);



%
% Set powers, phases, and frequencies
%
[P,ph,freq] = deal(zeros(numel(t),2));
for nn = 1:numPulses
    tc = t0 + (nn-1)*T;



    %
    % Set powers
    %
    if strcmpi(PulseType,'Gaussian') == 1
        P(:,1) = P(:,1) + power1(nn).*I_factor(:,1).*exp(-(t - tc).^2/fwhm.^2);
        P(:,2) = P(:,2) + power2(nn).*I_factor(:,1).*exp(-(t - tc).^2/fwhm.^2);
    elseif strcmpi(PulseType,'Square') == 1
        PulseEnd = find(abs(t-(tc + 1*width/2)) < dt/2,1,'first');
        PulseStart = find(abs(t-(tc - 1*width/2)) < dt/2,1,'last');
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
    % Set frequencies.
    % Need channel 1 frequency to be higher than channel
    % 2 frequency so that we use the lattice formed from retroreflecting
    % from the vibrationally isolated mirror
    %
    freq(idx,1) = DDSChannel.DEFAULT_FREQ + delta;
    freq(idx,2) = DDSChannel.DEFAULT_FREQ - delta;
end

freq(freq == 0) = DDSChannel.DEFAULT_FREQ;

%% Populate DDS values
for nn = 1:numel(dds)
    dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
end


end