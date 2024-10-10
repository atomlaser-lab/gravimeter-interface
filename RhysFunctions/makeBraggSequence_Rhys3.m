function makeBraggCompositePulse_Rhys3(dds,varargin)

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

AccelAmp = [];
NormalisedDisplacement = [];

w0 = 2.5e-3;
t0_effective = 0;
RampOnOff = 0;

% % % Pulse Fidelity
FreqError = 0;

% flag3
% % % Accelerometer Inputs
BW = 5e3;
ScaleFactor = 1;
Bias = 0;
White = 10;

if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 'pulsetype'
                PulseType = v;
            case 'freqerror'
                FreqError = v;
            case 'noisetype'
                NoiseType = v;
            case 'ramponoff'
                RampOnOff = v;
            case 'accelamp'
                AccelAmp = v;
            case 'normdisp'
                NormalisedDisplacement = v;
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


%% Select Pulse Type
% % % % Pulse Types
if strcmpi(PulseType,'Primitive')
    PulseAreas = [180];
    AllPhase = [0];
elseif strcmpi(PulseType,'Corpse')
    PulseAreas = [60 300 420];
    AllPhase = [0 180 0];
elseif strcmpi(PulseType,'BB1')
    PulseAreas = [180 360 180 180];
    AllPhase = [104.5 313.4 104.5 0];
elseif strcmpi(PulseType,'Knill')
    PulseAreas = [180 180 180 180 180];
    AllPhase = [240 210 300 210 240];
elseif strcmpi(PulseType,'Waltz')
    PulseAreas = [90 180 270];
    AllPhase = [0 180 0];
elseif strcmpi(PulseType,'N_90_360_90')
    PulseAreas = [90 360 90];
    AllPhase = [0 120 0];
elseif strcmpi(PulseType,'Scrofulous')
    PulseAreas = [180 180 180];
    AllPhase = [60 300 60];
elseif strcmpi(PulseType,'levitt')
    PulseAreas = [90 180 90];
    AllPhase = [90 0 90];
elseif strcmpi(PulseType,'N_90_240_90')
    PulseAreas = [90 240 90];
    AllPhase = [240 330 240];
else
    error('Not a listed pulse type')
end
% % % % 
% Convert Pulse Area into pulse durations, assuming a constant power
t_PiSegment = width;
t_segments = PulseAreas/180*t_PiSegment;
numSegments = numel(PulseAreas);
NumStd = 1.5;

% if isinteger(round(width/2*1e6,10)) == 0
%     warning('Pulse Duration/2 must be an integer')
%     return
% end


%% Creat Time Vector flag1
t = []; 
Phase = [];
for ii = 1:numSegments
    if ii == 1
        tPulseSegment = (-NumStd*t_segments(ii):dt:NumStd*t_segments(ii))';
        SegmentCentres(ii) = t0;
    else
        tPulseSegment = (-NumStd*t_segments(ii):dt:NumStd*t_segments(ii))';
        tPulseSegment = tPulseSegment - tPulseSegment(1) + t(end);
%         SegmentCentres(ii) = t0; %t0 + NumStd*t_segments(ii-1) + NumStd*t_segments(ii);
        SegmentCentres(ii) = SegmentCentres(ii - 1) + NumStd*t_segments(ii-1) + NumStd*t_segments(ii);
    end

    PhaseSegment = ones(size(tPulseSegment))*AllPhase(ii);
    Phase = [Phase PhaseSegment.'];

    fwhm(ii) = t_segments(ii)/(2*sqrt(log(2)));
    t = [t tPulseSegment.'];
    clear tPulseSegment PhaseSegment
end
t = t + t0;

%% Calculate intermediate values
recoil = const.hbar*k^2/(2*const.mRb*2*pi);
numPulses = numel(power);
% fwhm = width/(2*sqrt(log(2)));
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


%% Create Noise vectors
% % % % Create intensitity noise profile
if strcmpi(NoiseType,'acceleration')
    % % % Effective time vector
    dt_eff = 10e-6;
    t_end_eff = t(end) - (t0 - t0_effective);
    t_effective = (0:dt_eff:t_end_eff);

    t_effective_partial = t - (t0 - t0_effective);


    % % % Effective acceleration
    if isempty(AccelAmp) == 1 && isempty(NormalisedDisplacement) == 1
        AccelAmp = 0;
    elseif isempty(AccelAmp) == 0 && isempty(NormalisedDisplacement) == 1
        a_effective = AccelAmp*ones(size(t_effective));
    elseif isempty(AccelAmp) == 1 && isempty(NormalisedDisplacement) == 0
        r_final = NormalisedDisplacement*w0;
        t_final = (2*T + t0_effective);
        AccelAmp = 2*r_final/(t_final^2);
        a_effective = AccelAmp*ones(size(t_effective));
    else
        warning('Either two noise types or none at all')
    end

    if RampOnOff == 1
        % % % Measured acceleration
        SampleRate = 1/BW;
        SampleTimes = (0:SampleRate:t_end_eff);
        a_noBW = a_effective*ScaleFactor + Bias;

        % % % Add Bandwidth
        for ii = 1:length(SampleTimes)
            if SampleTimes(ii) == 0
                TempPos = 1;
            else
                TempPos = find(t_effective-SampleTimes(ii) < SampleRate/5);
            end
            MesPos(ii) = TempPos(end);
        end
        MesPos = unique(MesPos);


        for ii = 1:length(MesPos)
            if ii == 1
%                 a_Seen(MesPos(ii):MesPos(ii+1)) = mean(a_noBW(MesPos(ii):MesPos(ii+1)));
                a_BW(ii) = mean(a_noBW(MesPos(ii):MesPos(ii+1))) + normrnd(0,White);
                v_BW(ii) = 0;
                d_BW(ii) = 0;
            elseif ii == length(MesPos)
%                 a_Seen(MesPos(ii):length(t)) = mean(a_noBW(MesPos(ii-1):MesPos(ii)));
                a_BW(ii) = mean(a_noBW(MesPos(ii-1):MesPos(ii))) + normrnd(0,White);

                dt_Temp = (t_effective(MesPos(ii)) - t_effective(MesPos(ii-1)));
                v_BW(ii) = v_BW(ii-1) + a_BW(ii-1)*dt_Temp;
                d_BW(ii) = d_BW(ii-1) + v_BW(ii-1)*dt_Temp + 0.5*a_BW(ii-1)*dt_Temp^2;
            else
%                 a_Seen(MesPos(ii):MesPos(ii+1)) = mean(a_noBW(MesPos(ii):MesPos(ii+1)));
                a_BW(ii) = mean(a_noBW(MesPos(ii):MesPos(ii+1)))+ normrnd(0,White);

                dt_Temp = (t_effective(MesPos(ii)) - t_effective(MesPos(ii-1)));
                v_BW(ii) = v_BW(ii-1) + a_BW(ii-1)*dt_Temp;
                d_BW(ii) = d_BW(ii-1) + v_BW(ii-1)*dt_Temp + 0.5*a_BW(ii-1)*dt_Temp^2;
            end
        end
        % These are the quantaties seen at each sampletime (plot with
        % SampleTimes)

        % % % previous measurement is held until next measurement is made
% %         for ii = 1:length(t_effective)
% %             if ii == 1
% %                 v_Seen(ii) = 0;
% %                 d_Seen(ii) = 0;
% %             else
% %                 dt_temp = t_effective(ii) - t_effective(ii-1);
% %                 v_Seen(ii) = v_Seen(ii-1) + a_Seen(ii-1)*dt_temp;
% %                 d_Seen(ii) = d_Seen(ii-1) + v_Seen(ii-1)*dt_temp + 0.5*a_Seen(ii-1)*dt_temp^2;
% %             end
% %         end

        v_Seen = v_BW;
        d_Seen = d_BW;
        % These are the quantities seen at each time step during the "effective drop time"
        % plot against t_effective

        % % % Relate to true time vector
%         for ii = 1:length(t_effective_partial)
%             TempPos = find(t_effective-t_effective_partial(ii) < dt_eff/5);
%             RealPos(ii) = TempPos(end);
%         end
%         for ii = 1:length(t_effective_partial)
%             TempPos = find(abs(SampleTimes-t_effective_partial(ii)) < dt_eff/5);
%             RealPos(ii) = TempPos(end);
%         end
        for ii = 1:length(t_effective_partial)
%             if ii == 1
%                 RealPos(ii) = 1;
%             else
                TempPos = find((SampleTimes- t_effective_partial(ii)) < SampleRate);
                RealPos(ii) = TempPos(end);
%             end
        end
        d_Measured = d_Seen(RealPos).';
        % DDS instructions only exist when there's pulses. These are the
        % quantitites seen at those times
    end
    % % % Noise
    r = 0.5*AccelAmp*(t - (t0 - t0_effective)).^2;
    if RampOnOff == 1
        I_Noise_factor = exp(-2*(r - d_Measured).^2/w0^2);
    else
        I_Noise_factor = exp(-2*r.^2/w0^2);
    end

    % This is the effective displacement and the noise (including intenisty
    % ramp)
    % plot against t_effective_partial
%     figure(131);hold on
%     plot(t_effective_partial,I_Noise_factor)
5;

elseif strcmpi(NoiseType,'white')
    if AccelAmp > 1
        warning('Noise Amplitude cannot be greater than 100%. Noise set to 100%')
        AccelAmp = 1;
    end
    I_Noise_factor = AccelAmp*rand(1,length(t),1);
elseif strcmpi(NoiseType,'Offset')
    if AccelAmp < 0
        warning('Noise Amplitude less than zero means no pulses. Noise set to 0.1')
        AccelAmp = 0.1;
    end
    I_Noise_factor = AccelAmp;
elseif strcmpi(NoiseType,'white') == 0 && strcmpi(NoiseType,'acceleration')
    warning('Type must be "white" or "acceleration"')
    return
end


%
% Set powers, phases, and frequencies
%
[P,ph,freq] = deal(zeros(numel(t),2));
for nn = 1:numSegments
%     tc = t0 + (nn-1)*T + max((nn-2),0)*Tasym;
    tc = SegmentCentres(nn);
    idx = (t - t0) > (nn-1-0.5)*T;

    P(:,1) = P(:,1) + (I_Noise_factor.*(power1*exp(-(t - tc).^2/fwhm(nn).^2))).';
    P(:,2) = P(:,2) + (I_Noise_factor.*(power2*exp(-(t - tc).^2/fwhm(nn).^2))).';

%     ph(idx,2) = appliedPhase;

    freq(idx,1) = DDSChannel.DEFAULT_FREQ + mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order)*4*recoil) + (FreqError*recoil)/2e6;
    freq(idx,2) = DDSChannel.DEFAULT_FREQ - mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order)*4*recoil);
end


% ph()
figure(4312);clf
plot((t-t0)*1e6,P*1e3)
hold on
plot((t-t0)*1e6,Phase)

freq(freq == 0) = DDSChannel.DEFAULT_FREQ;
%  flag2
%% Populate DDS values
for nn = 1:numel(dds)
    dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
end


end