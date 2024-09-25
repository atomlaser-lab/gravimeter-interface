function makeBraggSequence_Rhys(dds,varargin)

%% Set up variables and parse inputs
f = 384.224e12;
k = 2*pi*f/const.c;
t0 = 10e-3;
width = 30e-6;
T = 1e-3;
Tasym = 0;
appliedPhase = 0;
power = 0.05*[1,2,1];
power1 = [];power2 = [];
chirp = 2*k*9.795/(2*pi);
order = 1;
start_order = 0;
mirrorSwitch = 1;
NoiseType = 'acceleration';
NoiseAmp = 0;
w0 = 2.5e-3;
t0_effective = 0;

if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 'noisetype'
                NoiseType = v;
            case 'noiseamp'
                NoiseAmp = v;
            case 'beamradius'
                w0 = v;
            case 't0_effective'
                t0_effective = v;
            case 't0'
                t0 = v;
            case 't'
                T = v;
            case 'dt'
                dt = v;
            case 'tasym'
                Tasym = v;
            case 'width'
                width = v;
            case {'appliedphase','phase'}
                appliedPhase = v;
            case 'power'
                power = v;
                if any(power < 0)
                    error('Power needs a value between 0 and 1.');
                elseif power > 1
                    error('Power needs a value between 0 and 1.');
                else
                    power = v;
                end
            case 'chirp'
                chirp = v;
            case 'f'
                f = v;
                k = 2*pi*f/const.c;
            case 'k'
                k = v;
            case 'power1'
                power1 = v;
            case 'power2'
                power2 = v;
            case 'order'
                if v == 0
                    error('Bragg order should be non-zero!');
                elseif round(v) ~= v
                    error('Bragg order must be an integer!');
                else
                    order = v;
                end
            case 'start_order'
                 if round(v) ~= v
                     error('Bragg start order must be an integer!');
                else
                    start_order = v;
                end
            case 'mirror'
                switch lower(v)
                    case 'isolated'
                        mirrorSwitch = 1;
                    case {'mot','rigid'}
                        mirrorSwitch = -1;
                    otherwise
                        error('''Mirror'' option can only be either ''isolated'' or ''mot''');
                end
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end

%% Conditions on the time step and the Bragg order
if width<50e-6
    dt=1e-6;
else
    intermediatewidth=width*10^6;
    dt = ceil(intermediatewidth/50)*10^-6;
end
     
%% Calculate intermediate values
recoil = const.hbar*k^2/(2*const.mRb*2*pi);
numPulses = numel(power);
fwhm = width/(2*sqrt(log(2)));
detuning=start_order*const.hbar*k^2/const.mRb;

if isempty(power1)
    power1 = power;
end
if isempty(power2)
    power2 = power;
end

if numel(appliedPhase) == 0
    tmp = zeros(1,numPulses);
    tmp(end) = appliedPhase;
    appliedPhase = tmp;
end


%% Create vectors
tPulse = (-5*width:dt:5*width)';
t = repmat(tPulse,1,numPulses);
for  nn = 1:numPulses
    t(:,nn) = t(:,nn) + t0 + (nn-1)*T + max((nn-2),0)*Tasym;
end
t = t(:);

% % % % Create intensitity noise profile
if strcmpi(NoiseType,'acceleration')
    r_final = NoiseAmp*w0;
    t_final = (2*T + t0_effective);
    a_effective = 2*r_final/(t_final^2);
%     a_effective = NoiseAmp;

    ScaleFactor = 1.05;
    Bias = 200e-6;
    a_measure = a_effective*ScaleFactor + Bias;
    a_ramp = a_measure;
%     a_ramp(1:numel(tPulse)) = a_ramp(1);
%     a_ramp(numel(tPulse)+1:2*numel(tPulse)) = a_ramp(numel(tPulse) + 1);
%     a_ramp(2*numel(tPulse)+1:3*numel(tPulse)) = a_ramp(numel(tPulse)*2 + 1);
    
    t_effective = t - (t0 - t0_effective);
    r = 0.5*a_effective*(t - (t0 - t0_effective)).^2;
    r_measure = 0*0.5*a_ramp*(t - (t0 - t0_effective)).^2;
    I_Noise_factor = exp(-2*(r - r_measure).^2/w0^2);
    5;
% %     r_final = NoiseAmp*w0;
% %     t_final = (2*T + t0_effective);
% %     a_effective = 2*r_final/(t_final^2);
% % %     a_effective = NoiseAmp;
% % 
% %     ScaleFactor = 1.05;
% %     Bias = 200e-6;
% %     Bandwidth = 5e3;
% % 
% %     t_effective = t - (t(1) - t0_effective);
% %     a_infBandwidth = a_effective*ScaleFactor + Bias;
% %     SampleTime = 0:1/Bandwidth:t_effective(end);

% % 
% %     differences = abs(t_effective - SampleTime);
% %     AllCases = differences < 1/(2*Bandwidth);
% %     ValidCases = find(any(AllCases,1));
% % 
% % %     SamplePos = find
% % 
% %     a_ramp = a_infBandwidth;
% % 
% % 
% % %     a_ramp(1:numel(tPulse)) = a_ramp(1);
% % %     a_ramp(numel(tPulse)+1:2*numel(tPulse)) = a_ramp(numel(tPulse) + 1);
% % %     a_ramp(2*numel(tPulse)+1:3*numel(tPulse)) = a_ramp(numel(tPulse)*2 + 1);
% %     
% %     t_effective = t - (t0 - t0_effective);
% %     r = 0.5*a_effective*(t - (t0 - t0_effective)).^2;
% %     r_measure = 0.5*a_ramp*(t - (t0 - t0_effective)).^2;
% %     I_Noise_factor = exp(-2*(r - r_measure).^2/w0^2);

elseif strcmpi(NoiseType,'white')
    if NoiseAmp > 1
        warning('Noise Amplitude cannot be greater than 100%. Noise set to 100%')
        NoiseAmp = 1;
    end
    I_Noise_factor = NoiseAmp*rand(1,length(t),1);
elseif strcmpi(NoiseType,'Offset')
    if NoiseAmp < 0
        warning('Noise Amplitude less than zero means no pulses. Noise set to 0.1')
        NoiseAmp = 0.1;
    end
    I_Noise_factor = NoiseAmp;
elseif strcmpi(NoiseType,'white') == 0 && strcmpi(NoiseType,'acceleration')
    warning('Type must be "white" or "acceleration"')
    return
end


%
% Set powers, phases, and frequencies
%
[P,ph,freq] = deal(zeros(numel(t),2));
for nn = 1:numPulses
    tc = t0 + (nn-1)*T + max((nn-2),0)*Tasym;
    idx = (t - t0) > (nn-1-0.5)*T;

    P(:,1) = P(:,1) + I_Noise_factor.*(power1(nn)*exp(-(t - tc).^2/fwhm.^2));
    P(:,2) = P(:,2) + I_Noise_factor.*(power2(nn)*exp(-(t - tc).^2/fwhm.^2));

    ph(idx,2) = appliedPhase(nn);

    freq(idx,1) = DDSChannel.DEFAULT_FREQ + mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order)*4*recoil);
    freq(idx,2) = DDSChannel.DEFAULT_FREQ - mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order)*4*recoil); 
end

freq(freq == 0) = DDSChannel.DEFAULT_FREQ;

%% Populate DDS values
for nn = 1:numel(dds)
    dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
end


end