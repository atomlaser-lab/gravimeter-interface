function sq = makeBEC_rhys(varargin)
% opt.LoadMagTrap_status = 0; opt.LoadOpticalTrap_status =0; opt.OpticalEvaporation_status = 0; opt.MagEvaporation_status = 0;
% opt.LoadMagTrap_status = 1; opt.LoadOpticalTrap_status =1; opt.OpticalEvaporation_status = 1; opt.MagEvaporation_status = 1;

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







%% Initialise Sequence: Set default values
sq = initSequence;

% % % The Raman DDS uses version 1.81 while Bragg DDS uses version 1.82
%  Insert new calibration data
% sq.dds(1).calibrationData.amp = ;
% sq.dds(1).calibraitonData.optical_power = ;


% % Version 1.81 requiers 
% sq.dds(1).rfscale = 3.1;
% sq.dds(2).rfscale = 2.15;
% Bragg Calibration data used in initSequence.
% For now I will simply load/set my own calibration data here. I'll create
% another object later
sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
calibData = load('RamanAOMData_formatted.mat');
sq.dds(1).calibrationData = calibData.data_ch1;
sq.dds(2).calibrationData = calibData.data_ch2;

Ch1_PMax = 59e-3; % @ 33 dBm
Ch2_PMax = 36.5e-3; % @ 33 dBm

I_ratio = 1;
Ch1_PDesire = 10e-3;

Ch2_PDesire = Ch1_PDesire*I_ratio;
Ch1_Pratio = Ch1_PDesire/Ch1_PMax;
Ch2_Pratio = Ch2_PDesire/Ch2_PMax;

if Ch1_PDesire > Ch1_PMax || Ch1_PDesire > Ch2_PMax
    error('You want more Raman power than you have')
end

%% Common Run options/functions
% functions
convert = RunConversions;
UD = @(x) x* -0.2031 + 2.707; %this converts the input value in amps to the appropriate voltage
DelayTime = @(X) 0.0107*X + 0.5776; % t in ms, this estimates the delay time needed between images for a camera pixel width of X
CentrePixelPos = @(t)2048 - (0.9138*t^2 + 0.1275*t + 670); %t in ms, this estimates vertical pixel that the centre of the cloud will be on in the image

% options
opt.misc.detuning2 = opt.detuning; % AOM cannot change fast enough
imageVoltage2 = convert.imaging(opt.misc.detuning2); % detuning for F = 1 atoms
imageVoltage = convert.imaging(opt.detuning); % detuning for F = 2 atoms

% misc options
RamanF1PumpTest = 0; % do you want to pump atoms into F = 1 after PGC?
RamanImagingSequence = 0; % do you want to use a 2-state imaging sequence?
Manifold = 1; %do you want to count atoms in the F=1 or F=2 state?

%% MOT
% % % Turn on dipoles at start of run
if opt.LoadOpticalTrap_status == 1
    sq.find('50w ttl').set(1);
    sq.find('25w ttl').set(1);
    sq.find('50w amp').set(convert.dipole50(26)); % 22
    sq.find('25w amp').set(convert.dipole25(19)); % 17
else
    sq.find('50w ttl').set(0);
    sq.find('25w ttl').set(0);
end
sq.find('Imaging Freq').set(convert.imaging(opt.detuning));
sq.find('Top Repump Shutter').set(0);


% % % Set 2D MOT Parameters
% Digital
sq.find('2D MOT Amp TTL').set(1);
sq.find('Push Amp TTL').set(1);
% Analogue
sq.find('Push Freq').set(9.5); % 9.5
sq.find('2D MOT Freq').set(7.75); %7.75

% % Set 3D Mot Parameters
% Digital
sq.find('3D MOT Amp TTL').set(1);
sq.find('Repump Amp TTL').set(1);
sq.find('MOT Coil TTL').set(1);
% Analogue
sq.find('3d coils').set(convert.mot_coil(1.6)); %1.6
sq.find('3D MOT Amp').set(5);
sq.find('3D MOT Freq').set(convert.mot_freq(-15)); %-15
sq.find('Repump Amp').set(9); %9
sq.find('Repump Freq').set(convert.repump_freq(-1.5)); %-1.5
sq.find('Liquid Crystal Repump').set(7);

% % % Mag Biases
sq.find('Bias E/W').set(0);
sq.find('Bias N/S').set(0);
sq.find('Bias U/D').set(UD(13)); %13, 0 V seems optimal (cancel ion pump field); however, check influence on PGC

% % % Duration
TMOT = opt.MOT_LoadTime;
sq.delay(TMOT);


%% Compressed MOT stage
if opt.CMOT_status == 1
    % % % Goal: produce a smaller/denser cloud

    % stop blowing hot atoms into the MOT
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);
    % Reduce re-radiation pressure (i.e. detune beams)
    sq.find('3D MOT freq').set(convert.mot_freq(-25)); % 26
    sq.find('repump freq').set(convert.repump_freq(-5)); %-8
    sq.find('3D MOT Amp').set(4.5); % 4.5
    sq.find('Repump Amp').set(9); %9
    % Increase/Decrease mag field to compress trapped atoms
    sq.find('3D coils').set(convert.mot_coil(1.3)); %1.45 % max OD at 1.7, max atom number at 1
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);
    sq.find('bias u/d').set(UD(13)); %13.4285
    % Duration
    Tcmot = 10e-3;
    sq.delay(Tcmot);
end

%% PGC stage
if opt.PGC_status == 1
    % Duration/Ramp
    Tpgc = 40*1e-3;
    t = linspace(0,Tpgc,100);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    % Detune/reduce light
    sq.find('3D MOT Amp').after(t,f(4.5,3.9)); %3.9
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-70))); %-71
    sq.find('Repump Amp').set(7); %7
    sq.find('repump freq').set(convert.repump_freq(-8.75)); %-9

    % Have B = 0 at the atoms
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),0));
    sq.find('bias u/d').set(UD(0));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);

    sq.delay(Tpgc);
    if RamanF1PumpTest == 0 && opt.LoadMagTrap_status == 0
        sq.find('repump amp ttl').set(0);
    end
end


if RamanF1PumpTest == 1
    % Pump atoms into the F = 1 state
    % i.e. turn repump off
    %     T = opt.params; %5*0.15e-3
    %     T = 1*0.15e-3; %5*0.15e-3
    if T == 0

    else
        sq.find('3d coils').set(convert.mot_coil(1.9));
        sq.find('repump amp ttl').set(0);
        sq.find('liquid crystal repump').set(-2.22);
        %     sq.find('3D MOT Freq').set(convert.mot_freq(-18));
        sq.find('3D MOT Amp').set(1);
        sq.delay(T);

        %     % Pump atoms into the F = 2 state
        %          % i.e. turn trapping off
        %     T = 30e-3;
        %     sq.find('3D MOT Amp TTL').set(0);
        %     sq.find('repump freq').set(convert.repump_freq(0));
        %     sq.find('Repump Amp').set(9);
        %     sq.delay(T);
    end
end

%% Load into magnetic trap
if opt.LoadMagTrap_status == 1
    % % % Pump atom into F = 1 manifold
    % Turn off repump
    Tdepump = 1e-3;
    sq.find('3D MOT Freq').set(convert.mot_freq(-70)); %-71
    sq.find('3D MOT Amp').set(4.5); %3.9
    sq.find('repump amp ttl').set(0);
    sq.find('Top repump shutter').set(1);
    sq.find('liquid crystal repump').set(-2.22);
    sq.delay(Tdepump);

    % % % Mag Load
    % Ramp mag field on
    tMLoad = 2*1e-3;
    t = linspace(0,tMLoad,100);
    f = @(vi,vf) sq.linramp(t,vi,vf);
    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').after(t,f(convert.mot_coil(5),convert.mot_coil(6.5)));

    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(2.5);
    sq.find('bias u/d').set(UD(0));

    % Turn light off
    sq.find('liquid crystal bragg').set(-3.9);
    sq.find('3D mot amp ttl').set(0);

    if opt.MagEvaporation_status == 0
        sq.delay(200e-3);
    end
end

%% Microwave evaporation
if opt.MagEvaporation_status ==1
    evapRate = 10; %10
    InitialFreq = 30; %30, 40
    FinalFreq = 4; % 12, 4
    Tevap = (InitialFreq-FinalFreq)/evapRate;

    t = linspace(0,Tevap,100);
    sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,InitialFreq,FinalFreq)));
    sq.find('mw amp ttl').set(1);
    sq.delay(Tevap);
end

%% Weaken trap while MW frequency fixed
if opt.LoadOpticalTrap_status == 1
    Trampcoils = 200*1e-3;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),convert.mot_coil(1))); % 1
    sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),0));
    sq.delay(Trampcoils);

    if opt.OpticalEvaporation_status == 0
        sq.find('bias e/w').set(0);
        sq.find('bias n/s').set(0);
        sq.find('bias u/d').set(0);
        sq.find('3d coils').set(0);
        sq.delay(100e-3);
    end
end

%% Optical evaporation
if opt.OpticalEvaporation_status == 1
    % Ramp down magnetic trap in 1 s

    Trampcoils = .8;
    t = linspace(0,Trampcoils,101);
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.find('mw amp ttl').anchor(sq.find('3d coils').last).before(100e-3,0);
    sq.find('mot coil ttl').at(sq.find('3d coils').last,0);
    % %
    % % At the same time, start optical evaporation
    % %
    sq.delay(30e-3);
    Tevap = 3;
    t = linspace(0,Tevap,200);
    FinalPower = 5.2;
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(FinalPower),0.5)); %4.8
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(FinalPower+0.25),0.5));

    sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),0)); %10 as end value?
    sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias u/d').values(end),0));
    sq.delay(Tevap);
end



%% Drop atoms
% sq.find('mot coil ttl').set(0);
% sq.find('3D Coils').set(convert.mot_coil(0));
% sq.delay(35e-3 - opt.tof);

timeAtDrop = sq.time; %Store the time when the atoms are dropped for later
sq.anchor(timeAtDrop);

sq.find('2D MOT Amp TTL').set(0);
sq.find('Push Amp TTL').set(0);

sq.find('3D mot amp ttl').set(0);
if opt.mw.enable(1) == 0 && opt.raman.OnOff == 0
    sq.find('bias e/w').before(50e-3,0); %200 ms before
end
sq.find('bias n/s').before(50e-3,0);
sq.find('bias u/d').before(50e-3,0);
sq.find('mw amp ttl').set(0);
sq.find('mot coil ttl').set(0);
sq.find('3D Coils').set(convert.mot_coil(0));
sq.find('25w ttl').set(0);
sq.find('50w ttl').set(0);
sq.find('50w amp').set(convert.dipole50(0));
sq.find('25w amp').set(convert.dipole25(0));

%% microwave transfer
if opt.mw.enable(1) == 1
    sq.anchor(timeAtDrop);

    % Inputs
        % microwave frequencies set out of run
        % no amplitude control
    MicrowaveDelay = 8e-3; %8
    MicrowaveDuration = 350*1e-6; %372,620
    
    % Set bias
    Delay = 500*1e-3;
    sq.find('bias e/w').before(Delay,10);
    sq.find('bias u/d').before(Delay,0); %%% This could be wrong. Previous operation was at zero VOLTS not amps
    sq.find('bias n/s').set(0);
    
    % Microwave Transfer
    sq.find('R&S list step trig').set(1);
    sq.delay(MicrowaveDelay); % delay to prevent state-changing collisions
    sq.find('state prep ttl').set(1);
    sq.delay(MicrowaveDuration);
    sq.find('state prep ttl').set(0);

    % return bias to zero
    sq.find('bias e/w').after(1e-3,0);
end



%% Interfeormetry

% % % if opt.raman.OnOff == 1
% % %     TimeBeforeDrop = 1.5;
% % % 
% % %     % Move the timing anchor
% % %     sq.anchor(timeAtDrop - TimeBeforeDrop);
% % % 
% % %     % timing
% % %     InterferometryDelay = 200e-3;        
% % %     PulseWidth = 100e-3;    
% % %     dt = 20e-3;
% % %     T_int = 200e-3;
% % %     
% % % 
% % %     % Trigger DDS 
% % %     TriggerDelay = 10e-3;
% % % %     sq.find('dds trig').before(TriggerDelay,1);
% % % %     sq.find('dds trig').after(TriggerDelay,0); %MOGLabs DDS triggers on falling edge
% % % %     sq.find('dds trig').after(TriggerDelay,1);
% % %     sq.find('Raman DDS Trig').before(TriggerDelay,1);
% % %     sq.find('Raman DDS Trig').after(TriggerDelay,0); %MOGLabs DDS triggers on falling edge
% % %     sq.find('Raman DDS Trig').after(TriggerDelay,1);
% % %     sq.ddsTrigDelay = timeAtDrop - TimeBeforeDrop;
% % % 
% % % 
% % %     % Inputs
% % %     Closer = 8e-3 - InterferometryDelay;    
% % %     chirp = 25.106258428e6;
% % %     Tasym =0;
% % %     k = 22.731334388721734;
% % % 
% % % %     delta = opt.params;
% % %     delta = -4.05;
% % %     
% % %     MakePulseSequence_Rhys_new(sq.dds,'k',k,'t0',InterferometryDelay,'T',T_int,'width',PulseWidth,'dt',dt,...
% % %         'phase',[0,0,0],'chirp',chirp,'delta',delta,...
% % %         'power1',1*[Ch1_Pratio,0,0],'power2',1*[Ch2_Pratio,0,0],'PulseType','Square','NumPulseWidths',1,'I_factor',1,'I_ratio',I_ratio);
% % % 
% % % end

if opt.raman.OnOff == 1
    % % % Inputs
    % Timing     
    Closer = 0e-3;
    InterferometryDelay = 8e-3 - Closer;
    T_int = 1e-3;
    PulseWidth = 500e-3;    
    dt = 50e-3;
    % Pulse Parameters
    chirp = 25.106258428e6;
    k = 22.731334388721734;
    delta = opt.params;
%     delta = 4;

    % Move the timing anchor
    sq.anchor(timeAtDrop);

    % Set bias
    Delay = 500*1e-3;
    sq.find('bias e/w').before(Delay,10);
    sq.find('bias u/d').before(Delay,0); %%% This could be wrong. Previous operation was at zero VOLTS not amps
    sq.find('bias n/s').set(0);

    % Trigger DDS
    TriggerDelay = 10e-3;
    sq.find('Raman DDS Trig').before(TriggerDelay,1);
    sq.find('Raman DDS Trig').after(TriggerDelay,0); %MOGLabs DDS triggers on falling edge
    sq.find('Raman DDS Trig').after(TriggerDelay,1);
    sq.ddsTrigDelay = timeAtDrop;
    

    % % % Add Intensity Noise
    I_factor = MakeIntensityNoise('type','acceleration','amp',0,'width',30e-6,'dt',1e-6,'NumPulseWidths',1);
    I_factor = 1;

    MakePulseSequence_Rhys_new(sq.dds,'k',k,'t0',InterferometryDelay,'T',T_int,'width',PulseWidth,'dt',dt,...
        'phase',[0,0,0],'chirp',chirp,'delta',delta,...
        'power1',1*[Ch1_Pratio,0,0],'power2',1*[Ch2_Pratio,0,0],'PulseType','Square','NumPulseWidths',1,'I_factor',1,'I_ratio',I_ratio);

    sq.find('bias e/w').after(Delay + InterferometryDelay + 1e-3,0);


    % % %     Test to see DDS is turning on
%         % Ch1 turn on for 1 second
%         sq.dds(1).set(110,1,0);
%         sq.dds(1).after(1,110,1,0);
%         sq.dds(1).after(1e-6,110,0,0);
%         % Ch 2 turn on for 1 second
%         sq.dds(2).set(110,0,0);
%         sq.dds(2).after(1,110,0,0);
%         sq.dds(2).after(1e-6,110,0,0);
end


%% Stern-Gerlach
if opt.mw.enable_sg == 1
    if opt.mw.enable(1) == 1
        Delay = MicrowaveDelay;
    elseif opt.raman.OnOff == 1
        Delay = InterferometryDelay + Closer;
    else
        Delay = 8e-3;
    end

    % inputs
    Tsg = 20e-3;    
    SternGerlachDelay = 2e-3;
    MaxCurrent = 4*1.5;
    
    sq.anchor(timeAtDrop + Delay + SternGerlachDelay);

    
    t = linspace(0,Tsg/2,20);
    
    sq.find('mot coil ttl').set(1);
    sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,MaxCurrent))); % ramp on
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0))); % ramp off
    sq.delay(Tsg);
    sq.find('mot coil ttl').set(0);
%     sq.find('3d coils').set(convert.mot_coil(0));
    sq.find('3d coils').set(0);

end





%% Imaging stage
%
% Image the atoms.  Reset the pointer for the whole sequence to when
% the atoms are dropped from the trap.  This means that the
% time-of-flight (tof) used in makeImagingSequence is now the delay
% from the time at which the atoms are dropped to when the first
% imaging pulse occurs
%
Abs_Analysis_parameters.camera = evalin('base', 'Abs_Analysis_parameters.camera');
sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 2;   %Set camera acquisition delay to be 2 s less than when image is taken
if strcmpi(Abs_Analysis_parameters.camera,'in-trap') || strcmpi(Abs_Analysis_parameters.camera,'drop 2')
    if RamanImagingSequence == 1
        makeRamanImagingSequence(sq,'type',Abs_Analysis_parameters.camera,...
            'tof1',opt.tof,'repump delay',0.5e-3,'tof2',opt.tof+opt.misc.tof2,'cycle time',150e-3,...
            'repump Time',200e-6,'pulse Delay',10e-6,'pulse time',10e-6,...
            'imaging freq',imageVoltage,'image freq2',imageVoltage2,'repump freq',4.3,'includeDarkImage',true);
    else
        makeImagingSequence(sq,'type',Abs_Analysis_parameters.camera,'tof',opt.tof,...
            'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
            'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.3,...
            'manifold',Manifold,'includeDarkImage',true,'cycle time',150e-3);
    end
elseif strcmpi(Abs_Analysis_parameters.camera,'drop 3') || strcmpi(Abs_Analysis_parameters.camera,'drop 4')
    makeFMISequence(sq,'tof',opt.tof,'offset',30e-3,'duration',100e-3,...
        'imaging freq',imageVoltage,'manifold',1);
end




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