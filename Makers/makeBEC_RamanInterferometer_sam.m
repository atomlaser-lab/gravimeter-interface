function sq = makeBEC_RamanInterferometer(varargin)
% Check if any variable is an instance of SequenceOptions
% Get all variable names in the workspace
allVarNames = evalin('base', 'who');

% Check if any variable is an instance of SequenceOptions
seqOptVarIndices = cellfun(@(varName) isa(evalin('base', varName), 'SequenceOptions'), allVarNames);
seqOptVarNames = allVarNames(seqOptVarIndices);

if isempty(seqOptVarNames)
    if  nargin ==1 && ~exist(varargin{1},'var')
        error('define the sequence option variable name first');
    end
    % Create a new instance of SequenceOptions if it does not exist
    opt = SequenceOptions;
    cprintf('Keywords','Initialising Sequence Options\n')
    assignin('base', 'opt', opt);

else
    if numel(seqOptVarNames) > 1
        error('Multiple SequenceOptions found in the workspace');
    end
    % Use the existing instance of SequenceOptions
    opt = evalin('base', seqOptVarNames{1});
end

% Check input arguments and update options accordingly
if nargin == 1
    if ~isa(varargin{1}, 'SequenceOptions')
        error('If using only one argument, it must be of type SequenceOptions');
    end
elseif nargin > 1
    if ~isa(varargin{1}, 'SequenceOptions')
        error('First argument must be of type SequenceOptions');
    end
    opt.replace(varargin{1});
    opt.set(varargin{2:end});
end

if nargout == 0
    % Search for the opt variable in the base workspace
    varExists = evalin('base', 'exist(''opt'', ''var'')');
    if varExists
        opt = evalin('base', 'opt');
    end
else
    varargout{1} = opt;
end

%% Create a conversion object to handle conversions to volts
convert = RunConversions;
imageVoltage = convert.imaging(opt.detuning);

%% Initialize sequence - defaults should be handled here
sq = initSequence;

% % Bragg Calibration data used in initSequence.
sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
calibData = load('RamanAOM_24062024');
sq.dds(1).calibrationData = calibData.data_ch1;
sq.dds(2).calibrationData = calibData.data_ch2;
% ch1 max = 33.68 dBm, ch2 max = 34.35 dBm

timeAtStart = sq.time;
sq.find('Imaging Freq').set(convert.imaging(opt.detuning));

if opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    sq.find('50w ttl').set(1);
    sq.find('25w ttl').set(1);
    sq.find('50w amp').set(convert.dipole50(35 - 7.75)); %22 %35
    sq.find('25w amp').set(convert.dipole25(20 - 3.5)); %17 %22 %20
end

%% Set up the MOT loading values
sq.find('liquid crystal repump').set(7);
sq.find('3d coils').set(RunConversions.mot_coil(1.35));
sq.find('3D MOT Freq').set(convert.mot_freq(-19));
sq.find('3D MOT Amp').set(RunConversions.mot_power(0.75));
sq.find('Repump freq').set(convert.repump_freq(-1.6));
sq.find('Repump Amp').set(RunConversions.repump_power(0.972));
sq.find('MOT coil TTL').set(1);
sq.find('bias u/d').set(RunConversions.UD_A_to_V(12.4));
sq.find('bias e/w').set(0); % flipped gives more atoms
sq.find('bias n/s').set(0);

sq.delay(opt.MOT_LoadTime);

%% Compressed MOT stage
if opt.CMOT_status == 1
    % %Turn off the 2D MOT and push beam 10 ms before the CMOT stage
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);

    %reduce re-radiation pressure, and tighten the trap
    sq.find('3D MOT freq').set(convert.mot_freq(-23));
    sq.find('3D MOT Amp').set(RunConversions.mot_power(0.883));
    sq.find('3D coils').set(RunConversions.mot_coil(1.77 + 0.2));

    sq.find('repump freq').set(convert.repump_freq(1.5));
    sq.find('Repump Amp').set(RunConversions.repump_power(0.9575));

    sq.find('Bias N/S').set(sq.find('Bias N/S').values(end));
    sq.find('Bias E/W').set(sq.find('Bias E/W').values(end));
    sq.find('Bias U/D').set(RunConversions.UD_A_to_V(8));

    sq.delay(13*1e-3);
end

%% PGC stage
if opt.PGC_status == 1
    Tpgc = 20*1e-3;
    t = linspace(0,Tpgc,50);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    sq.find('3D MOT Amp').after(t,f(sq.find('3D MOT Amp').values(end),RunConversions.mot_power(0.7)));
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-67.5)));
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),RunConversions.mot_coil(0)));

    sq.find('repump freq').set(convert.repump_freq(6.5));
    sq.find('Repump Amp').set(RunConversions.repump_power(1));

    sq.find('bias u/d').set(RunConversions.UD_A_to_V(8));
    sq.find('bias e/w').set(0.1).after(Tpgc,0);
    sq.find('bias n/s').set(0).after(Tpgc,0);
    sq.delay(Tpgc);
end

%% Load into magnetic trap
if opt.LoadMagTrap_status == 1 && opt.JustMOT ~= 1
    % Turn off the repump field for optical pumping - 1 ms
    Tdepump = 1e-3;
    sq.find('repump amp ttl').set(0);
    sq.find('Top repump shutter').before(2e-3,1);
    sq.delay(Tdepump);

    sq.find('liquid crystal bragg').set(-3);
    sq.find('3D MOT Amp').set(0);
    sq.find('3D mot amp ttl').set(0);
    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').set(RunConversions.mot_coil(10.05));


    % % % Return UD bias to cancel ion pump field
    t = linspace(0,100*1e-3,50);
    f = @(vi,vf) sq.linramp(t,vi,vf);
    sq.find('bias u/d').after(t,f(sq.find('bias u/d').values(end),RunConversions.UD_A_to_V(13)));

    %     sq.delay(0.1); %test lifetime in mag trap
end

%% Microwave evaporation
if opt.MagEvaporation_status == 1 && opt.LoadMagTrap_status == 1 && opt.JustMOT ~= 1
    sq.delay(22*1e-3);
    evapRate = 12.5;
    evapStart = 40;
    evapEnd = 8;
    Tevap = (evapStart-evapEnd)/evapRate;
    t = linspace(0,Tevap,100);
    sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,evapStart,evapEnd)));
    sq.find('mw amp ttl').set(1);
    sq.delay(Tevap);
end

%% Weaken trap while MW frequency fixed
if opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    Trampcoils = 80e-3;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),RunConversions.mot_coil(5.45)));
    sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
    sq.delay(Trampcoils);

    %Hold as atoms fall out of trap
    sq.delay(80e-3);
end

%% Optical evaporation
if opt.OpticalEvaporation_status == 1 && opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    % Ramp down magnetic trap in 1 s
    Trampcoils = 0.9 - 0.1 + 0.1;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.find('mw amp ttl').anchor(sq.find('3d coils').last).before(100e-3,0);
    sq.find('mot coil ttl').at(sq.find('3d coils').last,0);
    % %
    % % At the same time, start optical evaporation
    % %
    sq.delay(5*1e-3);
    Tevap = 3 - 0.5;
    t = linspace(0,Tevap,100);
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(1.5 + 0.8),0.42)); %+0.8, -0.3 for small
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(1.4 + 0.8),0.7)); % +0.9

    sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
    sq.delay(Tevap);
    time_at_evap_end = sq.time;
end


%% In Trap Raman Transfer: |1,-1> -> |2,0>
% % % bean
% P1_max = 3.99; %mW %ch1 max at 34.14 dBm
% P2_max = 29.15; % mW %ch2 max at 34.35 dBm
% BlowDuration = 1.5*1e-3;
% 
% % % Time
% InTrapDelay1 = 100e-3;
% dt = 1e-6;
% RamanPulseWidth = 10*1e-6;
% TriggerDelay = 10e-6;
% 
% % % Power
% P2onP1 = 12;
% P_total = 30;
% 
% % % Detuning
% delta = -20;
% 
% % % Bias fields
% BiasEW = 0; %2
% BiasUD = 4; % 4
% BiasNS = 0;
% RamanBiasRampTime = 50e-3;
% RamanBiasDelay = 60e-3; %100
% 
% % % % Intermediate values
% % convert power to AOM setting
% P1 = P_total/(1+P2onP1);
% P2 = (P_total*P2onP1)/(1+P2onP1);
% Ch1_AOMSetting = P1/P1_max;
% Ch2_AOMSetting = P2/P2_max;
% % account for DDS error
% RamanPulseWidth = RamanPulseWidth - dt;
% 
% % % % Make Pulse
% if opt.StatePrep == 1
%     % Bias
%     sq.anchor(time_at_evap_end - InTrapDelay1 - RamanBiasDelay);
%     t_on = linspace(0,RamanBiasRampTime,100);
%     t_off = linspace(0,10e-3,11);
%     sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD));
%     sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS));
%     sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW));
%     sq.anchor(time_at_evap_end - InTrapDelay1 + RamanPulseWidth + 100e-6);
%     BiasEndPoint = time_at_evap_end - InTrapDelay1 + RamanPulseWidth + 100e-6;
%     if opt.raman == 0
%         sq.find('Bias U/D').after(t_off,sq.minjerk(t_off,sq.find('bias U/D').values(end),0));
%         sq.find('Bias N/S').after(t_off,sq.minjerk(t_off,sq.find('bias N/S').values(end),0));
%         sq.find('Bias E/W').after(t_off,sq.minjerk(t_off,sq.find('bias E/W').values(end),0));
%     end
% 
%     % Trigger DDS
%     sq.anchor(time_at_evap_end - InTrapDelay1 - TriggerDelay);
%     sq.find('DDS Trig').set(0).after(10e-3,1);
%     sq.ddsTrigDelay = time_at_evap_end - InTrapDelay1 - TriggerDelay;
% 
%     MakePulseSequence_Rhys(sq.dds,...
%         't0',TriggerDelay,'width',RamanPulseWidth,'dt',dt,...
%         'phase',[0 0 0],'delta',delta,...,
%         'power1',[Ch1_AOMSetting,0,0],...
%         'power2',[Ch2_AOMSetting,0,0],...
%         'PulseShape','Square');
% 
%     if opt.mw.analyze(1) ~= 1
%         % % % Turn on repump for in-trap blow away of remaining F = 1 atoms
%         sq.anchor(time_at_evap_end);
%         sq.find('Repump Amp TTL').before(InTrapDelay1 - RamanPulseWidth,1).after(BlowDuration,0);
%         sq.find('Top Repump Shutter').before(InTrapDelay1 + 3e-3 - RamanPulseWidth,0).after(BlowDuration + 1e-3,1);
%         sq.find('repump freq').before(InTrapDelay1 - RamanPulseWidth,convert.repump_freq(0));
%         sq.find('Repump Amp').before(InTrapDelay1 - RamanPulseWidth,10);
%     end
%     sq.anchor(time_at_evap_end);
% end
%% In Trap MW Transfer: |1,-1> -> |2,0> -> |1,0>
% % % Inputs
ExtraDelay = 0*1e-3;

% pulse 1
MW1Duration = 750*1e-6;

InTrapDelay1 = 100e-3 + ExtraDelay;
BlowDuration1 = 1.5*1e-3;

BiasEW = 4; %4 when single power supply
BiasUD = RunConversions.UD_A_to_V(13);
BiasNS = 0;
MWBiasRampTime = 50e-3;
MWBiasDelay = 100e-3; %100


% pulse 2
MW2Duration = 210*1e-6;
InTrapDelay2 = 50e-3 + ExtraDelay;%10e-3
BlowDuration2 = 12*1e-6;

if opt.mw.enable(1) == 1 % % % Transfer |1,-1> -> |2,0>    % Bias
    sq.anchor(time_at_evap_end - InTrapDelay1 - MWBiasDelay);
    t_on = linspace(0,MWBiasRampTime,100);
    t_off = linspace(0,30e-3,20);
%     sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD));
    sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS));
    sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW));
    sq.anchor(time_at_evap_end - InTrapDelay1 + MW1Duration + 100e-6);
    BiasEndPoint = time_at_evap_end - InTrapDelay1 + MW1Duration + 100e-6;
    if opt.raman == 0
        sq.delay(1e-3);
%         sq.find('Bias U/D').after(t_off,sq.minjerk(t_off,sq.find('bias U/D').values(end),0));
        sq.find('Bias N/S').after(t_off,sq.minjerk(t_off,sq.find('bias N/S').values(end),0));
        sq.find('Bias E/W').after(t_off,sq.minjerk(t_off,sq.find('bias E/W').values(end),0));
    end




    t_off = linspace(0,40e-3,40);
    sq.anchor(time_at_evap_end);
    sq.find('state prep ttl').before(InTrapDelay1,1).after(MW1Duration,0);
    if opt.mw.analyze(1) ~= 1
        % % % Turn on repump for in-trap blow away of remaining F = 1 atoms
        sq.find('Repump Amp TTL').before(InTrapDelay1 - MW1Duration,1).after(BlowDuration1,0);
        sq.find('Top Repump Shutter').before(InTrapDelay1 + 3e-3 - MW1Duration,0).after(BlowDuration1 + 1e-3,1);
        sq.find('repump freq').before(InTrapDelay1 - MW1Duration,convert.repump_freq(0));
        sq.find('Repump Amp').before(InTrapDelay1 - MW1Duration,10);
    end
end
if opt.mw.enable(2) == 1 && opt.mw.analyze(1) ~= 1 % % % Transfer |2,0> -> |1,0>
    sq.anchor(time_at_evap_end);
    sq.find('R&S list step trig').before(InTrapDelay2+40e-3,0);
    sq.find('state prep ttl').before(InTrapDelay2,1).after(MW2Duration,0);
    if opt.mw.analyze(2) ~= 1
        % % % Turn on trap for in-trap blow away of remaining F = 2 atoms
        sq.find('3D MOT Amp TTL').before(InTrapDelay2 - MW2Duration,1).after(BlowDuration2,0);
        sq.find('3D MOT Amp').before(InTrapDelay2 - MW2Duration,5).after(BlowDuration2,0);
        sq.find('3D MOT Freq').before(InTrapDelay2 - MW2Duration,RunConversions.mot_freq(0)).after(BlowDuration2,RunConversions.mot_freq(-75));
    end
end

%% Drop atoms
timeAtDrop = sq.time;
sq.anchor(timeAtDrop);

sq.find('3D mot amp ttl').set(0);
if opt.raman == 0
    sq.find('bias e/w').before(200e-3,0);
    sq.find('bias n/s').before(200e-3,0);
%     sq.find('bias u/d').before(200e-3,0);
end
sq.find('mw amp ttl').set(0);
sq.find('mot coil ttl').set(0);
sq.find('3D Coils').set(convert.mot_coil(0));
sq.find('25w ttl').set(0);
sq.find('50w ttl').set(0);
sq.find('50w amp').set(convert.dipole50(0));
sq.find('25w amp').set(convert.dipole25(0));
sq.find('2D MOT Amp TTL').before(10e-3,0);
sq.find('push amp ttl').before(10e-3,0);
sq.find('repump amp ttl').before(10e-3,0);

%% Raman Stuff
% Flag
P1_max = 3.18; %mW %ch1 max at 34.14 dBm
P2_max = 37.15; % mW %ch2 max at 34.35 dBm
SixShots = 0;
CompositePulseOnOff = 0;
PulseType = PulseTypes.Waltz;

% % % % Inputs
Pulse1OnOff = 1;
Pulse2OnOff = 0;
Pulse3OnOff = 0;

% % Time
dt = 1e-6;
T_Sep = 0.1e-3;
RamanPulseWidth = 10*1e-6;
RamanPulseWidth2 = 2*RamanPulseWidth;
RamanTOF = 10e-6 + 0*1e-3;

% % Power
P2onP1 = 7;
P_total = 20.77;

% % Detuning
delta = -20 + 0*opt.params*1e-3;
k = 2*pi/const.c*384.229441689483e12;
chirp = 0*2*k*9.795/(2*pi); %Chirp in Hz/s.

% % Bias fields
BiasEW = 10;
BiasUD = RunConversions.UD_A_to_V(4.5);
BiasNS = 0;
RamanBiasRampTime = 50e-3;
RamanBiasDelay2 = 100e-3;


% % Phase
phi_1 = 0;
phi_2 = 0;
phi_3 = 137.74;

% % Intensity Ramp
Gravity = 9.8;
w0 = 2.65/2*1e-3*0.8;
RampOnOff = 0;



% % % % % % % % % Calculate power in each pulse & AOM setting
if opt.raman == 1
    if RampOnOff == 1
        Ramp_1 = exp((0.5*Gravity*(RamanTOF)^2)^2/w0^2);
        Ramp_2 = exp((0.5*Gravity*(RamanTOF + T_Sep)^2)^2/w0^2);
        Ramp_3 = exp((0.5*Gravity*(RamanTOF + 2*T_Sep)^2)^2/w0^2);
    else
        Ramp_1 = 1;
        Ramp_2 = 1;
        Ramp_3 = 1;
    end

    P1_pulse1 = Pulse1OnOff*Ramp_1*P_total/(1+P2onP1);
    P2_pulse1 = Pulse1OnOff*(Ramp_1*P_total*P2onP1)/(1+P2onP1);
    P1_pulse2 = Pulse2OnOff*Ramp_2*P_total/(1+P2onP1);
    P2_pulse2 = Pulse2OnOff*(Ramp_2*P_total*P2onP1)/(1+P2onP1);
    P1_pulse3 = Pulse3OnOff*Ramp_3*P_total/(1+P2onP1);
    P2_pulse3 = Pulse3OnOff*(Ramp_3*P_total*P2onP1)/(1+P2onP1);

    if P1_pulse1> P1_max || P1_pulse2 > P1_max || P1_pulse3 > P1_max
        error('Channel 1 has insufficient power')
    end
    if P2_pulse1 > P2_max || P2_pulse2 > P2_max || P2_pulse3 > P2_max
        error('Channel 2 has insufficient power')
    end
    Ch1_AOMSetting_pulse1 = P1_pulse1/P1_max;
    Ch2_AOMSetting_pulse1 = P2_pulse1/P2_max;
    Ch1_AOMSetting_pulse2 = P1_pulse2/P1_max;
    Ch2_AOMSetting_pulse2 = P2_pulse2/P2_max;
    Ch1_AOMSetting_pulse3 = P1_pulse3/P1_max;
    Ch2_AOMSetting_pulse3 = P2_pulse3/P2_max;
end

% % % % % % % % Correct for pulse timing error
% The dds adds an extra instruction equal to dt
% e.g. a 100 us pulse made with dt = 10 us is 110 us long
% e.g. a 100 us pulse made with dt = 100 us is 200 us long
% e.g. a 0 us pulse made with dt = 1 is 1 us long
% Hence:
RamanPulseWidth = RamanPulseWidth - dt;
RamanTOF = RamanTOF - dt;
T_Sep = T_Sep + dt;
% % % % % % % %



% % % % Sequence

% Stern-Gerlach pulse to test MW transfer
if opt.mw.enable_sg == 1
    %
    % Apply a Stern-Gerlach pulse to separate states based on magnetic
    % moment.  A ramp is used to ensure that the magnetic states
    % adiabatically follow the magnetic field
    %
    % % %     % works with 25 ms tof
    SGDelay = 3e-3; %3
    SGDuration = 8e-3; %8

    sq.anchor(timeAtDrop + SGDelay);
    sq.find('mot coil ttl').set(1);
    t = linspace(0,SGDuration,40);
    sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,5)));%5
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.delay(2*SGDuration);
    sq.find('mot coil ttl').set(0);
    sq.find('3d coils').set(convert.mot_coil(0));
end

% Raman transfer from |1,0> to |2,0>
if opt.raman == 1

    NumPulses = Pulse1OnOff + Pulse2OnOff + Pulse3OnOff;

    RamanPulseStartTime = timeAtDrop + RamanTOF;
    RamanPulseEndTime = RamanPulseStartTime + (NumPulses - 1)*T_Sep + RamanPulseWidth*(Pulse1OnOff + Pulse3OnOff) + RamanPulseWidth2*Pulse2OnOff;
    RamanBiasStartTime = RamanPulseStartTime - RamanBiasDelay2 - RamanBiasRampTime;

    if RamanBiasDelay2 >= InTrapDelay1
        sq.anchor(BiasEndPoint);
        t_on = linspace(0,(RamanPulseStartTime - BiasEndPoint)*0.8,100);
        sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD));
        sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS));
        sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW));
    else
        sq.anchor(RamanBiasStartTime);
        t_on = linspace(0,RamanBiasRampTime,100);
        t_off = linspace(0,10e-3,11);
        sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD));
        sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS));
        sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW));
    end

    sq.anchor(RamanPulseEndTime + 100e-6);
    sq.find('Bias U/D').after(t_off,sq.minjerk(t_off,sq.find('bias U/D').values(end),0));
    sq.find('Bias N/S').after(t_off,sq.minjerk(t_off,sq.find('bias N/S').values(end),0));
    sq.find('Bias E/W').after(t_off,sq.minjerk(t_off,sq.find('bias E/W').values(end),0));

    % % % Phase lock
    sq.anchor(timeAtDrop + RamanTOF - 10e-3);
    sq.find('Raman Phase Lock').set(1).after(10e-3 + 4*RamanPulseWidth + 2*T_Sep + 1e-3,0);



    % % % % Make Raman pulse(s)
    sq.anchor(timeAtDrop);
    if opt.StatePrep == 0
        sq.find('DDS Trig').set(0).after(10e-3,1);
        sq.ddsTrigDelay = timeAtDrop;
    end

%     if CompositePulseOnOff == 1
%         MakeCompositePulse_Rhys(sq.dds,...
%             'PulseType',PulseType,'delta',delta,...
%             'P1_max',P1_max,'P2_max',P2_max,'P_pi',P_total,'P_rat',P2onP1,...
%             't0',RamanTOF,'width',RamanPulseWidth,'dt',dt);
%     else
%         MakePulseSequence_Rhys(sq.dds,...
%             't0',RamanTOF,'T',T_Sep,'width',RamanPulseWidth,'width2',RamanPulseWidth2,'dt',dt,...
%             'phase',[phi_1,phi_2,phi_3],'delta',delta,...,
%             'power1',[Ch1_AOMSetting_pulse1,Ch1_AOMSetting_pulse2,Ch1_AOMSetting_pulse3],...
%             'power2',[Ch2_AOMSetting_pulse1,Ch2_AOMSetting_pulse2,Ch2_AOMSetting_pulse3],...
%             'PulseShape','Square','chirp',chirp);
%     end
    
    T = 0.05e-3; %interferometer spacing
    numPulses = 3; %pulse count
    % following arrays must contain at least the pulse number of points
    pulsewidth = 1e-6*[4,8,4];
    power1 = 0.5*[1,2,1];
    power2 = 0.5*[1,2,1];
    appliedPhase = [0,0,180];

    %frequnecy

    delta = 20; %Mhz
    freq1 = 110 + delta/4;
    freq2 = 110 - delta/4;
    
    %time array at pulse points, 2 per pulse
    t = zeros(1,numPulses*2);
    for i =1:numPulses
        t(2*i-1) = (i-1)*T; %start of pulse is at nT 
        t(2*i) = (i-1)*T+pulsewidth(i); %end of pulse after width
    end
    
    %interleve powers with zeros
    B =zeros(size(power1));
    C = [power1;B];
    power1 = C(:)';
    B =zeros(size(power2));
    C = [power2;B];
    power2 = C(:)';

    %duplicate phases
    appliedPhase = repelem(appliedPhase,2);
    
    % Set powers, phases, and frequencies
    %
    [P,ph,freq] = deal(zeros(numel(t),2));
    %
    % Set powers
    P(:,2) = power2;
    P(:,1) = power1;
    %
    % Set phases
    %
   ph(:,2) = appliedPhase;


    %
    % Set frequencies
    %
    freq(:,1) = freq(:,1) +freq1;
    freq(:,2) = freq(:,2) +freq2;   
    
    % Populate DDS values
    for nn = 1:numel(sq.dds)
        if nn == 1
            sq.dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));        
        elseif nn == 2
            sq.dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
        end
end

end

%% Imaging stage
%
% Image the atoms.  Reset the pointer for the whole sequence to when
% the atoms are dropped from the trap.  This means that the
% time-of-flight (tof) used in makeImagingSequence is now the delay
% from the time at which the atoms are dropped to when the first
% imaging pulse occurs
%
% Abs_Analysis_parameters.camera = evalin('base', 'Abs_Analysis_parameters.camera');
sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 2;   %Set camera acquisition delay to be 2 s less than when image is taken
% if strcmpi(Abs_Analysis_parameters.camera,'in-trap') || strcmpi(Abs_Analysis_parameters.camera,'drop 2')
if opt.TwoStateImaging == 1
    makeRamanImagingSequence(sq,'type','in-trap',...
        'tof1',opt.tof,'repump delay',0.5e-3,'tof2',opt.tof+opt.misc.tof2,'cycle time',150e-3,...
        'repump Time',200e-6,'pulse Delay',10e-6,'pulse time',10e-6,...
        'imaging freq',imageVoltage,'image freq2',imageVoltage,'repump freq',4.3,'includeDarkImage',true);
else
    makeImagingSequence(sq,'type','in-trap','tof',opt.tof,...
        'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
        'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.2,...
        'manifold',1,'includeDarkImage',true,'cycle time',150e-3);
end

if SixShots == 1
    sq.anchor(sq.latest);
    sq.delay(150e-3);
    sq.dds(1).set(DDSChannel.DEFAULT_FREQ,0,0);
    sq.find('cam trig').set(1).after(10e-6,0);
    sq.delay(10e-6);
    sq.dds(1).set(DDSChannel.DEFAULT_FREQ,1e-3/P1_max,0);
    sq.delay(150e-6);
    sq.dds(1).set(DDSChannel.DEFAULT_FREQ,0,0);
    sq.delay(150e-3);
    sq.dds(2).set(DDSChannel.DEFAULT_FREQ,0,0);
    sq.find('cam trig').set(1).after(10e-6,0);
    sq.delay(10e-6);
    sq.dds(2).set(DDSChannel.DEFAULT_FREQ,1.5e-3/P2_max,0);

    sq.delay(150e-6);
    sq.dds(2).set(DDSChannel.DEFAULT_FREQ,0,0);
end

sq.find('Liquid Crystal Repump').at(timeAtDrop,7);

%% Automatic start
%If no output argument is requested, then compile and run the above
%sequence
if nargout == 0
    r = RemoteControl;
    r.upload(sq.compile);
    r.run;
else
    varargout{1} = sq;
end
end