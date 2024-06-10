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


% sq.find('Cam Trig').set(1).after(1e-3,0);
sq.find('Cam Trig').set(0);
sq.find('DDS Trig').set(0);

% % Bragg Calibration data used in initSequence.
% % For now I will simply load/set my own calibration data here. I'll create
% % another object later
sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
% calibData = load('raman-aom-calibration-29-04-2024.mat');
calibData = load('BraggDDS_RamanAOM');
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

sq.find('3D MOT Freq').set(convert.mot_freq(-14));
sq.find('3D MOT Amp').set(5);

sq.find('Repump freq').set(convert.repump_freq(-1.1));
sq.find('Repump Amp').set(8);

sq.find('MOT coil TTL').set(1);
sq.find('bias u/d').set(0);
sq.find('bias e/w').set(0);
sq.find('bias n/s').set(0);
sq.delay(opt.MOT_LoadTime);

%% Compressed MOT stage
if opt.CMOT_status == 1
    % %Turn off the 2D MOT and push beam 10 ms before the CMOT stage
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);

    %reduce re-radiation pressure, and tighten the trap
    sq.find('3D MOT freq').set(convert.mot_freq(-22 -1));
    sq.find('3D MOT Amp').set(4.4);
    sq.find('3D coils').set(RunConversions.mot_coil(1.47 + 0.3));

    %     sq.find('3D MOT freq').set(convert.mot_freq(-14 -5));
    %     sq.find('3D MOT Amp').set(RunConversions.mot_power(0.884));
    %     sq.find('3D coils').set(RunConversions.mot_coil(1.35 + 0.8));

    sq.find('repump freq').set(convert.repump_freq(-8 - 1));
    sq.find('Repump Amp').set(7.5);

    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);
    sq.find('bias u/d').set(0);

    sq.delay(13*1e-3);
end

%% PGC stage
if opt.PGC_status == 1
    Tpgc = 20*1e-3;
    t = linspace(0,Tpgc,50);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    %Smooth ramps for these parameters
    sq.find('3D MOT Amp').after(t,f(sq.find('3D MOT Amp').values(end),3.8));
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-69  +1 + 0.5)));
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),RunConversions.mot_coil(0)));

    sq.find('repump freq').set(convert.repump_freq(-9 + 0.1));
    sq.find('Repump Amp').set(6.25);

    sq.find('bias u/d').set(3 - 2);
    sq.find('bias e/w').set(0).after(Tpgc,0); %"Flipped" EW bias
    sq.find('bias n/s').set(10).after(Tpgc,0);
    %     sq.find('bias n/s').before(10e-3,0).after(Tpgc,0);


    sq.delay(Tpgc);

    % Turn off the repump field for optical pumping - 1 ms
    Tdepump = 1e-3;
    sq.find('repump amp ttl').set(0);
    sq.find('Top repump shutter').before(2e-3,1);
    sq.delay(Tdepump);
end

%% Load into magnetic trap
if opt.LoadMagTrap_status == 1 && opt.JustMOT ~= 1
    sq.find('liquid crystal bragg').set(-3);
    sq.find('3D MOT Amp').set(0);
    sq.find('3D mot amp ttl').set(0);
    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').set(RunConversions.mot_coil(10.05));
    % sq.delay(3.0); %test lifetime in mag trap
end

%% Microwave evaporation
if opt.MagEvaporation_status == 1 && opt.LoadMagTrap_status == 1 && opt.JustMOT ~= 1
    sq.delay(22*1e-3);
    evapRate = 14;
    evapStart = 40;
    evapEnd = 9.5;
    Tevap = (evapStart-evapEnd)/evapRate;
    t = linspace(0,Tevap,100);
    sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,evapStart,evapEnd)));
    sq.find('mw amp ttl').set(1);
    sq.delay(Tevap);
end

%% Weaken trap while MW frequency fixed
if opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    Trampcoils = 80e-3; %180
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),RunConversions.mot_coil(4.45)));
    sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),0));
    sq.delay(Trampcoils);

    %Hold as atoms fall out of trap
    sq.delay(80e-3);
end

%% Optical evaporation
if opt.OpticalEvaporation_status == 1 && opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    % Ramp down magnetic trap in 1 s
    Trampcoils = 0.9 - 0.1;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.find('mw amp ttl').anchor(sq.find('3d coils').last).before(100e-3,0);
    sq.find('mot coil ttl').at(sq.find('3d coils').last,0);
    % %
    % % At the same time, start optical evaporation
    % %
    sq.delay(5*1e-3);
    Tevap = 2.5 + 0.5;
    t = linspace(0,Tevap,100);
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(1.5 + 0.4),0.42));
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(1.4 + 0.4),0.7));

    sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias u/d').values(end),0));
    sq.delay(Tevap);
    time_at_evap_end = sq.time;
end



%% In Trap MW Transfer: |1,-1> -> |2,0> -> |1,0>
% % % Inputs
% pulse 1

% MW1Duration = 3*550*1e-6;
MW1Duration = 525*1e-6;
InTrapDelay1 = 100e-3;
BlowDuration1 = 1.5*1e-3;
% pulse 2
MW2Duration = 200*1e-6;
InTrapDelay2 = 50e-3;%10e-3
BlowDuration2 = 15*1e-6;

if opt.mw.enable(1) == 1 % % % Transfer |1,-1> -> |2,0>
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
if opt.mw.enable(2) == 1 % % % Transfer |2,0> -> |1,0>
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
sq.find('bias e/w').before(200e-3,0);
sq.find('bias n/s').before(200e-3,0);
sq.find('bias u/d').before(200e-3,0);
sq.find('mw amp ttl').set(0);
sq.find('mot coil ttl').set(0);
sq.find('3D Coils').set(convert.mot_coil(0));
sq.find('25w ttl').set(0);
sq.find('50w ttl').set(0);
sq.find('50w amp').set(convert.dipole50(0));
sq.find('25w amp').set(convert.dipole25(0));


sq.find('raman dds trig').before(20e-3,0).after(1e-3,1);
sq.find('dds trig').before(20e-3,0).after(1e-3,1);
sq.ddsTrigDelay = sq.time - 20e-3;
%% Raman Stuff
% Flag (easy for ctrl f)
% % % % Inputs
Pulse1OnOff = 1;
Pulse2OnOff = 0;

phi_1 = 0;
phi_2 = 0;
T_Sep = 5e-3;

% Bias fields
BiasEW = 10; %9.5
BiasUD = 0; %2
BiasNS = 9; % 0
RamanBiasDelay = 40*1e-3;

RamanTOF = 0*1e-6;
RamanPulseWidth = 50*1e-6; % 9
% RamanPulseWidth = opt.params*1e-6; % 9
dt = 1e-6;

% delta = -20 + opt.params*1e-3;
delta = -20 - 0*1e-3;

P2onP1 = (7/1);
% P_total = 4;
P_total = opt.params;

P1_max = 4.155; %mW
P2_max = 25.47; % mW

Sideband_badPol_max = 5.5e-6; %uW
Carrier_badPol_max = 66.6e-6; %uW


P1 = P_total/(1+P2onP1);
P2 = (P_total*P2onP1)/(1+P2onP1);

Ch1_AOMSetting = P1/P1_max;
Ch2_AOMSetting = P2/P2_max;
% Ch1_AOMSetting = 1;
% Ch2_AOMSetting = 1;

if P1 > P1_max
    error('Channel 1 has insufficient power')
end
if P2 > P2_max
    error('Channel 2 has insufficient power')
end





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

    triggerDelay = 1e-3;
    TriggerDuration = 10e-3;

    % % % Pulse 1
    sq.anchor(timeAtDrop + RamanTOF);
    sq.find('DDS Trig').before(TriggerDuration + triggerDelay,1);
    sq.find('DDS Trig').after(TriggerDuration,0);
    sq.ddsTrigDelay = timeAtDrop + RamanTOF - triggerDelay;

    sq.anchor(timeAtDrop + RamanTOF - RamanBiasDelay);
    t_on = linspace(0,RamanBiasDelay,30);
    t_off = linspace(0,5e-3,15);
    sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD)).after(t_off,sq.minjerk(t_off,sq.find('bias U/D').values(end),0));
    sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS)).after(t_off,sq.minjerk(t_off,sq.find('bias N/S').values(end),0));
    sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW)).after(t_off,sq.minjerk(t_off,sq.find('bias E/W').values(end),0));

    % % % % Make Raman pulse(s)
    sq.anchor(timeAtDrop);
    if Pulse1OnOff == 1 && Pulse2OnOff == 0
        MakePulseSequence_Rhys(sq.dds,'t0',RamanTOF,'T',T_Sep,'width',RamanPulseWidth,'dt',dt,...
            'phase',[phi_1,phi_2,0],'delta',delta,...
            'power1',[Ch1_AOMSetting,0,0],'power2',[Ch2_AOMSetting,0,0],'PulseType','Square');
    elseif Pulse1OnOff == 0 && Pulse2OnOff == 1
        MakePulseSequence_Rhys(sq.dds,'t0',RamanTOF,'T',T_Sep,'width',RamanPulseWidth,'dt',dt,...
            'phase',[phi_1,phi_2,0],'delta',delta,...
            'power1',[0,Ch1_AOMSetting,0],'power2',[0,Ch2_AOMSetting,0],'PulseType','Square');
    elseif Pulse1OnOff == 1 && Pulse2OnOff == 1
        MakePulseSequence_Rhys(sq.dds,'t0',RamanTOF,'T',T_Sep,'width',RamanPulseWidth,'dt',dt,...
            'phase',[phi_1,phi_2,0],'delta',delta,'rampAmp',9.8,...,
            'power1',[Ch1_AOMSetting,Ch1_AOMSetting,0],'power2',[Ch2_AOMSetting,Ch2_AOMSetting,0],'PulseType','Square');
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
% elseif strcmpi(Abs_Analysis_parameters.camera,'drop 3') || strcmpi(Abs_Analysis_parameters.camera,'drop 4')
%     makeFMISequence(sq,'tof',opt.tof,'offset',30e-3,'duration',100e-3,...
%         'imaging freq',imageVoltage,'manifold',1);
% end
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