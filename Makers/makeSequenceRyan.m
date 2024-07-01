function sq = makeSequenceRyan(varargin)
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
sq.find('DDS Trig').set(1);
sq.find('Raman DDS Trig').set(1);

% % Bragg Calibration data used in initSequence.
% % For now I will simply load/set my own calibration data here. I'll create
% % another object later
% sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
% sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
% calibData = load('BraggDDS_RamanAOM');
calibData = load('RamanAOMHex_RamanDDS_14062024');
sq.dds(1).calibrationData = calibData.data_ch1;
sq.dds(2).calibrationData = calibData.data_ch2;
sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;

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
sq.find('3D MOT Amp').set(4.2);
sq.find('3D MOT Freq').set(convert.mot_freq(-17.25 + 3.9 - 0.5));
sq.find('Repump freq').set(convert.repump_freq(-1.5 + 0.5 - 0.1));
sq.find('MOT coil TTL').set(1);
sq.find('3d coils').set(0.15);
sq.find('bias u/d').set(0);
sq.find('bias e/w').set(5);
sq.find('bias n/s').set(0);
sq.delay(opt.MOT_LoadTime);

%% Compressed MOT stage
if opt.CMOT_status == 1
    % %Turn off the 2D MOT and push beam 10 ms before the CMOT stage
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);

    %Increase the cooling and repump detunings to reduce re-radiation
    %pressure, and weaken the trap
    sq.find('3D MOT freq').set(convert.mot_freq(-22 -1));
    sq.find('repump freq').set(convert.repump_freq(-8 - 1));
    sq.find('3D coils').set(RunConversions.mot_coil(1.47 + 0.3));
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
    sq.find('3D MOT Amp').after(t,f(4.2,3.8 - 0.2));
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-69  +1))); %66 for 12s
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),RunConversions.mot_coil(0)));
    sq.find('repump freq').set(convert.repump_freq(-9 + 0.3));

    sq.find('bias u/d').set(3 - 2);
    %     sq.find('bias e/w').set(10).after(Tpgc,0); %"normal" EW bias
    sq.find('bias e/w').set(0).after(Tpgc,0); %"Flipped" EW bias
    sq.find('bias n/s').set(0).after(Tpgc,0);

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
    evapRate = 14 + 2 + 3;
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
    Trampcoils = 180e-3;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),RunConversions.mot_coil(3.45 + 1)));
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
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(1.5 + 0.8),0.42));
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(1.5 + 0.8),0.7));

    sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias u/d').values(end),0));
    sq.delay(Tevap);
    time_at_evap_end = sq.time;
end



%% In Trap MW Transfer: |1,-1> -> |2,0> -> |1,0>
% % % Inputs
% pulse 1
MW1Duration = 580*1e-6;
InTrapDelay1 = 100e-3;
BlowDuration1 = 1.5*1e-3;
% pulse 2
MW2Duration = 200*1e-6;
InTrapDelay2 = 50e-3;%10e-3
BlowDuration2 = 10*1e-6;

if opt.mw.enable(1) == 1 % % % Transfer |1,-1> -> |2,0>
    sq.anchor(time_at_evap_end - InTrapDelay1);
    sq.find('state prep ttl').set(1);
    sq.delay(MW1Duration);
    sq.find('state prep ttl').set(0);
    if ~opt.mw.analyze(1)
        % % % Turn on repump for in-trap blow away of remaining F = 1 atoms
        sq.find('Repump Amp TTL').set(1);
        sq.find('repump freq').set(convert.repump_freq(0));
        sq.find('Repump Amp').set(10);
        sq.find('Top Repump Shutter').before(3e-3,0);
        sq.delay(BlowDuration1);
        sq.find('repump amp ttl').set(0);
        sq.find('repump amp').set(0);
        sq.find('Top Repump Shutter').at(sq.time,1);
    end
end
if opt.mw.enable(2) == 1 % % % Transfer |2,0> -> |1,0>
    sq.anchor(time_at_evap_end - InTrapDelay2);
    sq.find('R&S list step trig').before(40e-3,0);    
    sq.find('state prep ttl').set(1);
    sq.delay(MW2Duration);
    sq.find('state prep ttl').set(0); % % % This was set to 1
    sq.find('R&S list step trig').at(sq.time,1);
    if ~opt.mw.analyze(2)
        % % % Turn on trap for in-trap blow away of remaining F = 2 atoms
        sq.find('3D MOT Amp TTL').set(1);
        sq.find('3D MOT Amp').set(5);
        sq.find('3D MOT Freq').set(convert.mot_freq(0));
        sq.delay(BlowDuration2);
        sq.find('3D MOT Amp TTL').set(0);
        sq.find('3D MOT Amp').set(0);
        sq.find('3D MOT Freq').set(convert.mot_freq(-75));
    end
end
sq.anchor(time_at_evap_end);

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
% % % Inputs
Raman = 0;

% RamanTOF = opt.params(3);
% RamanPulseWidth = opt.params(1)*1e-6;
RamanTOF = 2000e-6;
RamanPulseWidth = 150*1e-6; %pi = 300 for mag 3
dt = 1e-6;

P2onP1 = 7;
P_total = 5;

P1_max = 4.36; %mW
P2_max = 27.1; % mW

P1 = P_total/(1+P2onP1);
P2 = (P_total*P2onP1)/(1+P2onP1);

Ch1_AOMSetting = P1/P1_max;
Ch2_AOMSetting = P2/P2_max;
% Ch1_AOMSetting = 0.5;
% Ch2_AOMSetting = 0.1;

% delta = -20 + opt.params*1e-3;
delta = -20 - 0.756e-3;


% Bias fields
BiasUD = 1; %1
BiasNS = 10; % 10
BiasEW = 10; %10
RamanBiasDelay = 40*1e-3;


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

% Raman transfer from |2,0> to |1,0>
if Raman == 1
    sq.anchor(timeAtDrop + RamanTOF - RamanBiasDelay);
    t_on = linspace(0,RamanBiasDelay,30);
    t_off = linspace(0,5e-3,51);
%     sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD)).after(t_off,sq.minjerk(t_off,sq.find('bias U/D').values(end),0));
%     sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS)).after(t_off,sq.minjerk(t_off,sq.find('bias N/S').values(end),0));
%     sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW)).after(t_off,sq.minjerk(t_off,sq.find('bias E/W').values(end),0));

    sq.anchor(timeAtDrop + RamanTOF);
    t = (-4*dt):dt:(RamanPulseWidth + 4*dt);
    P1 = Ch1_AOMSetting*(t <= RamanPulseWidth & t >= 0);
    P2 = Ch2_AOMSetting*(t <= RamanPulseWidth & t >= 0);
    df = delta*ones(size(t));
    f1 = DDSChannel.DEFAULT_FREQ + df/4;
    f2 = DDSChannel.DEFAULT_FREQ - df/4;
    sq.dds(1).after(t,f1,P1,zeros(size(t)));
    sq.dds(2).after(t,f2,P2,zeros(size(t)));

%     sq.delay(RamanPulseWidth + 100e-6);
% 
%     sq.find('3D MOT Amp TTL').set(1);
%     sq.find('3D MOT Amp').set(5);
%     sq.find('3D MOT Freq').set(convert.mot_freq(0));
%     sq.delay(10e-6);
%     sq.find('3D MOT Amp TTL').set(0);
%     sq.find('3D MOT Amp').set(0);
%     sq.find('3D MOT Freq').set(convert.mot_freq(-75));
%     sq.delay(500e-6);
% 
%     t = (-4*dt):dt:(RamanPulseWidth2 + 4*dt);
%     P1 = Ch1Power*(t <= RamanPulseWidth2 & t >= 0);
%     P2 = Ch2Power*(t <= RamanPulseWidth2 & t >= 0);
%     df = delta2*ones(size(t));
%     f1 = DDSChannel.DEFAULT_FREQ + df/4;
%     f2 = DDSChannel.DEFAULT_FREQ - df/4;
%     sq.dds(1).after(t,f1,P1,zeros(size(t)));
%     sq.dds(2).after(t,f2,P2,zeros(size(t)));
   
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
        makeRamanImagingSequence(sq,'type','in trap',...
            'tof1',opt.tof,'repump delay',0.5e-3,'tof2',opt.tof+opt.misc.tof2,'cycle time',150e-3,...
            'repump Time',200e-6,'pulse Delay',10e-6,'pulse time',10e-6,...
            'imaging freq',imageVoltage,'image freq2',imageVoltage,'repump freq',4.3,'includeDarkImage',true);
    else
        makeImagingSequence(sq,'type','in trap','tof',opt.tof,...
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