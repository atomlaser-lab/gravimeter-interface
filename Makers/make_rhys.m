function sq = make_rhys(varargin)
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
% sq.find('Raman DDS Trig').set(1); %% DDS triggers on falling edge, so start on

% Bragg Calibration data used in initSequence.
% For now I will simply load/set my own calibration data here. I'll create
% another object later
sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
calibData = load('RamanAOMData_formatted.mat');
sq.dds(1).calibrationData = calibData.data_ch1;
sq.dds(2).calibrationData = calibData.data_ch2;

Ch1_PMax = 64e-3; % @ 33 dBm, Half waveplate @ 8 deg, Amp @
Ch2_PMax = 64e-3; % @ 33 dBm, Half waveplate @ 303 deg, Amp @

I_ratio = 1;
Ch1_PDesire = 10e-3;
Ch1_PDesire = Ch1_PMax;

Ch2_PDesire = Ch1_PDesire*I_ratio;
Ch1_Pratio = Ch1_PDesire/Ch1_PMax;
Ch2_Pratio = Ch2_PDesire/Ch2_PMax;
if Ch1_PDesire > Ch1_PMax || Ch1_PDesire > Ch2_PMax
    error('You want more Raman power than you have')
end

%% Common Run options/functions

% functions
convert = RunConversions;
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
sq.find('Imaging Freq').set(convert.imaging(opt.detuning));
% Start with dipoles on if loading into dipoles
if opt.LoadOpticalTrap_status == 1
    sq.find('50w ttl').set(1);
    sq.find('25w ttl').set(1);
    sq.find('50w amp').set(convert.dipole50(24.5)); % 22
    sq.find('25w amp').set(convert.dipole25(22)); % 19
else
    sq.find('50w ttl').set(0);
    sq.find('25w ttl').set(0);
end

% Goal: Collect as many atoms as possible
% % % % 2D MOT Light
sq.find('2D MOT Amp TTL').set(1);
sq.find('Push Amp TTL').set(1);

sq.find('2D MOT Freq').set(7.8);
sq.find('Push Freq').set(9.8 - 0.6);

% % % % 3D MOT Light
sq.find('3D MOT Amp TTL').set(1);
sq.find('Repump Amp TTL').set(1);
sq.find('MOT Coil TTL').set(1);
% % Trapping
sq.find('3D MOT Freq').set(convert.mot_freq(-18));
sq.find('3D MOT Amp').set(5);
sq.find('Liquid crystal Bragg').set(3);
% % Repump
sq.find('Repump Freq').set(convert.repump_freq(-1.5));
sq.find('Repump Amp').set(9);
sq.find('Liquid Crystal Repump').set(7);
%
% % % % Mag Field
sq.find('Bias E/W').set(0);
sq.find('Bias N/S').set(0);
sq.find('Bias U/D').set(0);
sq.find('3D Coils').set(convert.mot_coil(1.5));
% sq.find('3D Coils').set(convert.mot_coil(1.5 - 0.3));

% Duration
sq.delay(opt.MOT_LoadTime);

%% CMOT
if opt.CMOT_status == 1
    % Goal: Increase confinement with minimal loss
    % % % Stop blowing hot atoms in
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);

    % % % Reduce Re-radiation Pressure
    sq.find('3D MOT freq').set(convert.mot_freq(-18));
    sq.find('3D MOT Amp').set(5);

    sq.find('repump freq').set(convert.repump_freq(-8));
    sq.find('Repump Amp').set(7); %9

    % % % Increase Mag Field
    sq.find('3D Coils').set(convert.mot_coil(convert.mot_coil_reverse(sq.find('3D Coils').values(end)) - 0.5)); % test this with mag load
    sq.find('Bias E/W').set(0);
    sq.find('Bias N/S').set(0);
    sq.find('Bias U/D').set(0);

    % % %Duration
    sq.delay(15*1e-3);
end


%% PGC
if opt.PGC_status == 1
    % Duration/Ramp
    Tpgc = 35*1e-3;
    t = linspace(0,Tpgc,100);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    % Detune/reduce light
    sq.find('3D MOT Amp').after(t,f(sq.find('3D MOT Amp').values(end),3.6)); %3.9
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-71))); %-71
    sq.find('Repump Amp').set(5); %7
    sq.find('repump freq').set(convert.repump_freq(-9)); %-9

    % Have B = 0 at the atoms
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),0));
    sq.find('bias u/d').set(3); %4
    sq.find('bias e/w').set(0); %optimum in negatives
    sq.find('bias n/s').set(0); %optimum in negatives

    sq.delay(Tpgc);
    if RamanF1PumpTest == 0 && opt.LoadMagTrap_status == 0
        sq.find('repump amp ttl').set(0);
    end
end

%% Mag Load
if opt.LoadMagTrap_status == 1
    % % % Pump into |F=1, mf = -1>
    % Repump off
    sq.find('Repump Amp TTL').set(0);
    sq.find('Top repump shutter').set(1);
    sq.find('liquid crystal repump').set(-2.22);
    % Trapping amplitude
    sq.find('3D MOT Freq').set(convert.mot_freq(-71));
    sq.find('3D MOT Amp').set(4.5);
    sq.delay(0.6*1e-3);

    % % % Ramp on Mag Field
    t = linspace(0,2*1e-3,100);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').after(t,f(convert.mot_coil(3.5 + 1),convert.mot_coil(3.5 + 4 + 1)));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);
    sq.find('bias u/d').set(0);

    % % % Turn Trapping light Light Off
    sq.find('3D MOT Amp TTL').set(0);
    sq.find('Liquid crystal Bragg').set(-3.9);

    if opt.MagEvaporation_status == 0
        sq.delay(200*1e-3);
    end
end



%% Mag Evap
if opt.MagEvaporation_status == 1
    InitialFreq = 36; %33
    FinalFreq = 2;
    Rate = 11; % 10
    t_evap = (InitialFreq - FinalFreq)/Rate;
    t = linspace(0,t_evap,100);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    sq.find('mw amp ttl').set(1);
    sq.find('mw freq').after(t,f(convert.microwave(InitialFreq),convert.microwave(FinalFreq)));
    sq.delay(t_evap);
    sq.find('mw amp ttl').set(0);
end

%% Dipole Load
if opt.LoadOpticalTrap_status == 1
    % Loosen mag trap
    t_loosen = 50*1e-3;
    t = linspace(0,t_loosen,100);
    f = @(vi,vf) sq.minjerk(t,vi,vf);
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),convert.mot_coil(1 + 2.5)));
    sq.find('bias e/w').after(t,f(sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t,f(sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t,f(sq.find('bias u/d').values(end),0));
    sq.delay(t_loosen);
    % hold in dipoles
    if opt.OpticalEvaporation_status == 0
        sq.find('MOT Coil TTL').set(0);
        sq.find('3D coils').set(0);
        sq.delay(50*1e-3);
    end
end

%% Optical Evap
if opt.OpticalEvaporation_status == 1
    % Ramp down magnetic trap
    Trampcoils = 500*1e-3;
    t = linspace(0,Trampcoils,101);
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.find('mw amp ttl').anchor(sq.find('3d coils').last).before(100e-3,0);
    sq.find('mot coil ttl').at(sq.find('3d coils').last,0);
    % %
    % % At the same time, start optical evaporation
    % %
    sq.delay(50*1e-3);
    Tevap = 2.5;
    t = linspace(0,Tevap,200);
    FinalPower = 1.6; %1.6
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(FinalPower),0.7));
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(FinalPower),0.7));

    sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),0)); %10 as end value?
    sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias u/d').values(end),0));
    sq.delay(Tevap);
end

%% Drop
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

%% MW Transfer between mf states
if opt.mw.enable(1) == 1
    sq.anchor(timeAtDrop);

    % Inputs
    % microwave frequencies set out of run
    % no amplitude control
    MicrowaveDelay = 0e-3; %8
    MicrowaveDuration = 650*1e-6;
    MagDelay = 500*1e-3;

    % Set bias
    sq.find('bias e/w').before(MagDelay,10);
    sq.find('bias u/d').before(MagDelay,0); %%% This could be wrong. Previous operation was at zero VOLTS not amps
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
%% Raman Alignment
RamanAlignment = 1;
if RamanAlignment == 1 && opt.raman.OnOff ~= 1
    Ch2_Pratio = 0;
    Ch1_Pratio = 1;

    % % % Inputs
    % Timing
    TriggerDuration = 10e-3; % minimum of 10 ms needed
    triggerDelay = 1e-3; % minimum delay of 1 ms required
    Closer = 8e-3;
    InterferometryDelay = 8e-3 - Closer;
    T_int = 1e-3;

    TimeBeforeDrop = 20e-3;    %i.e. do pulse in trap
    PulseWidth = TimeBeforeDrop;
    dt = 1e-3;
    % Pulse Parameters
    chirp = 25.106258428e6;
    k = 22.731334388721734;
    delta = 5;

    % Set bias
    sq.anchor(timeAtDrop);
    MagDelay = 500*1e-3;
    sq.find('bias e/w').before(MagDelay,10);
    sq.find('bias u/d').before(MagDelay,0);
    sq.find('bias n/s').set(0);

    % Trigger DDS
    if mod(triggerDelay,1e-6) < 1e-6 && mod(triggerDelay,1e-6) ~= 0
        error('DDS Error: Trigger Delay requires DDS timing resolution less than 1 us')
    end
    sq.anchor(timeAtDrop + InterferometryDelay - triggerDelay);
    sq.find('Raman DDS Trig').before(TriggerDuration,1);
    sq.find('Raman DDS Trig').after(TriggerDuration,0);
    sq.ddsTrigDelay = timeAtDrop + InterferometryDelay - triggerDelay - TimeBeforeDrop;

    sq.anchor(timeAtDrop);
    MakePulseSequence_Rhys(sq.dds,'k',k,'t0',InterferometryDelay - TimeBeforeDrop,'T',T_int,'width',PulseWidth,'dt',dt,...
        'phase',[0,0,0],'chirp',chirp,'delta',delta,...
        'power1',1*[Ch1_Pratio,0,0],'power2',1*[Ch2_Pratio,0,0],'PulseType','Square');
end



%% Interferometry
% if opt.raman.OnOff == 1
% %     Ch2_Pratio = 0.1;
% %     Ch1_Pratio = 0.1;
% %         Ch2_Pratio = 0.;
% %         Ch1_Pratio = 0.;
% Ch1_Pratio = opt.params;
% Ch2_Pratio = opt.params;
% delta = 4.750;
% % % % Need to re-calibrate channels
% % % % Ch1_Pratio = Ch1_Pratio - 0.005; 
% % % % Ch2_Pratio = Ch2_Pratio - 0.015;
% % This will give 10 mW in each channel
% 
% 
%     % % % Inputs
%     % Timing
%     TriggerDuration = 10e-3; % minimum of 10 ms needed
%     triggerDelay = 1e-3; % minimum delay of 1 ms required
%     Closer = 8e-3;
%     InterferometryDelay = 8e-3 - Closer;
%     T_int = 1e-3;
%     PulseWidth = 10e-6;
%     dt = 1e-6;
%     % Pulse Parameters
%     chirp = 25.106258428e6;
%     k = 22.731334388721734;
%     delta = opt.params;
% 
%     % Set bias
%     sq.anchor(timeAtDrop);
%     MagDelay = 500*1e-3;
%     sq.find('bias e/w').before(MagDelay,10);
%     sq.find('bias u/d').before(MagDelay,0);
%     sq.find('bias n/s').set(0);
% 
%     % Trigger DDS
%     if mod(triggerDelay,1e-6) < 1e-6 && mod(triggerDelay,1e-6) ~= 0
%         error('DDS Error: Trigger Delay requires DDS timing resolution less than 1 us')
%     end
%     sq.anchor(timeAtDrop + InterferometryDelay - triggerDelay);
%     sq.find('Raman DDS Trig').before(TriggerDuration,1);
%     sq.find('Raman DDS Trig').after(TriggerDuration,0);
%     sq.ddsTrigDelay = timeAtDrop + InterferometryDelay - triggerDelay;
% 
%     sq.anchor(timeAtDrop);
%     MakePulseSequence_Rhys(sq.dds,'k',k,'t0',InterferometryDelay,'T',T_int,'width',PulseWidth,'dt',dt,...
%         'phase',[0,0,0],'chirp',chirp,'delta',delta,...
%         'power1',1*[Ch1_Pratio,0,0],'power2',1*[Ch2_Pratio,0,0],'PulseType','Square');
% 
%     sq.find('bias e/w').after(MagDelay + InterferometryDelay + 1e-3,0);
% 
% end


%% Stern-Gerlach
if opt.mw.enable_sg == 1
    if opt.mw.enable(1) == 1
        MagDelay = MicrowaveDelay;
    elseif opt.raman.OnOff == 1
        MagDelay = InterferometryDelay + Closer;
    else
        MagDelay = 8e-3;
    end

    % inputs
    Tsg = 20e-3;
    SternGerlachDelay = 2e-3;
    MaxCurrent = 3*1.5;

    sq.anchor(timeAtDrop + MagDelay + SternGerlachDelay);
    t = linspace(0,Tsg/2,20);
    sq.find('mot coil ttl').set(1);
    sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,MaxCurrent))); % ramp on
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0))); % ramp off
    sq.delay(Tsg);

    sq.find('mot coil ttl').set(0);
    sq.find('3d coils').set(0);
end


%% image
Abs_Analysis_parameters.camera = evalin('base', 'Abs_Analysis_parameters.camera');
% Abs_Analysis_parameters.camera = opt.imaging_type;
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