function sq = makeBEC_Bragg_NewCloud_backup(varargin)
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

sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_HEX_INTERP;
calibData = load('BraggPower_20082024');
sq.dds(1).calibrationData = calibData.data_ch1;
sq.dds(2).calibrationData = calibData.data_ch2;

timeAtStart = sq.time;
sq.find('Imaging Freq').set(convert.imaging(opt.detuning));

if opt.LoadOpticalTrap_status == 1
    sq.find('50w ttl').set(1);
    sq.find('25w ttl').set(1);
    sq.find('50w amp').set(convert.dipole50(16));
    sq.find('25w amp').set(convert.dipole25(15));
end

%% MOT
sq.find('liquid crystal repump').set(7);
sq.find('Top repump shutter').set(0);

sq.find('3D MOT Freq').set(convert.mot_freq(-18.5));
sq.find('3D MOT Amp').set(RunConversions.mot_power(1));
sq.find('Repump freq').set(convert.repump_freq(-1.5));
sq.find('Repump Amp').set(RunConversions.repump_power(1));

sq.find('MOT coil TTL').set(1);
sq.find('3d coils').set(RunConversions.mot_coil(1.35));
sq.find('bias u/d').set(convert.UD_A_to_V(13.8));
sq.find('bias e/w').set(1);
sq.find('bias n/s').set(0);

sq.delay(opt.MOT_LoadTime);

%% Compressed MOT stage
if opt.CMOT_status == 1
    PushDelay = 2*1e-3;
    sq.find('2D MOT Amp TTL').before(PushDelay,0);
    sq.find('push amp ttl').before(PushDelay,0);

    sq.find('3D MOT freq').set(convert.mot_freq(-23));
    sq.find('3D MOT Amp').set(RunConversions.mot_power(0.8));
    sq.find('repump freq').set(convert.repump_freq(-9));
    sq.find('Repump Amp').set(RunConversions.repump_power(1));
    
    sq.find('3D coils').set(RunConversions.mot_coil(1.5));
    sq.find('bias e/w').set(1);
    sq.find('bias n/s').set(0.5);
    sq.find('bias u/d').set(convert.UD_A_to_V(11)); %13

    Tcmot = 15*1e-3;
    sq.delay(Tcmot);
end

%% PGC stage
if opt.PGC_status == 1
    Tpgc = 20e-3;
    t = linspace(0,Tpgc,50);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    sq.find('3D MOT Amp').after(t,f(RunConversions.mot_power(1),RunConversions.mot_power(0.8)));
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-65)));
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),RunConversions.mot_coil(0.6)));
    sq.find('repump freq').set(convert.repump_freq(-9));

    sq.find('bias u/d').set(convert.UD_A_to_V(11));
    sq.find('bias e/w').set(1);
    sq.find('bias n/s').set(0.5);

    sq.delay(Tpgc);
end

%% Load into magnetic trap
if opt.LoadMagTrap_status == 1 && opt.JustMOT ~= 1
    Tdepump = 2e-3;
    sq.find('Top repump shutter').before(1.5e-3,1);
    sq.find('repump amp ttl').set(0);
    sq.find('liquid crystal repump').set(-2.522);
    sq.delay(Tdepump);

    sq.find('liquid crystal bragg').set(-3.9);
    sq.find('3D MOT Amp').set(0);
    sq.find('3D mot amp ttl').set(0);
    sq.find('3D coils').set(RunConversions.mot_coil(10));
end

%% Microwave evaporation
if opt.MagEvaporation_status == 1 && opt.LoadMagTrap_status == 1 && opt.JustMOT ~= 1
    sq.delay(22*1e-3);
    evapRate = 12;
    evapStart = 34;
    evapEnd = 8;
    Tevap = (evapStart-evapEnd)/evapRate;
    t = linspace(0,Tevap,100);
    sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,evapStart,evapEnd)));
    sq.find('mw amp ttl').set(1);
    sq.delay(Tevap);
end

%% Weaken trap while MW frequency fixed
if opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    Trampcoils = 50*1e-3;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),RunConversions.mot_coil(4.6)));
    sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),convert.UD_A_to_V(13.8)));
    sq.delay(Trampcoils);
end

%% Optical evaporation
if opt.OpticalEvaporation_status == 1 && opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    Trampcoils = 0.8;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),0));
    sq.find('mw amp ttl').anchor(sq.find('3d coils').last).before(100e-3,0);
    sq.find('mot coil ttl').at(sq.find('3d coils').last,0);

    sq.delay(10*1e-3);
    Tevap = 2.5;
    t = linspace(0,Tevap,100);
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(1.3 + 0.4),0.56));
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(1.3 + 0.4),0.72)); 
    sq.delay(Tevap);
    time_at_evap_end = sq.time;    
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


%% In Trap MW Transfer: |1,-1> -> |2,0> -> |1,0>
% % % % Inputs
% ExtraDelay = 0*1e-3;
% % pulse 1
% MW1Duration = 700*1e-6;
% 
% InTrapDelay1 = 100e-3 + ExtraDelay;
% BlowDuration1 = 1.5*1e-3;
% 
% BiasEW = 4;
% BiasUD = RunConversions.UD_A_to_V(13);
% BiasNS = 0;
% MWBiasRampTime = 50e-3;
% MWBiasDelay = 100e-3; %100
% 
% 
% % pulse 2
% MW2Duration = 210*1e-6;
% InTrapDelay2 = 50e-3 + ExtraDelay;%10e-3
% BlowDuration2 = 12*1e-6;
% 
% if opt.mw.enable(1) == 1 % % % Transfer |1,-1> -> |2,0>    % Bias
%     sq.anchor(time_at_evap_end - InTrapDelay1 - MWBiasDelay);
%     t_on = linspace(0,MWBiasRampTime,100);
%     t_off = linspace(0,30e-3,20);
%     %     sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD));
%     sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS));
%     sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW));
%     sq.anchor(time_at_evap_end - InTrapDelay1 + MW1Duration + 100e-6);
%     BiasEndPoint = time_at_evap_end - InTrapDelay1 + MW1Duration + 100e-6;
%     if opt.raman == 0
%         sq.delay(1e-3);
%         %         sq.find('Bias U/D').after(t_off,sq.minjerk(t_off,sq.find('bias U/D').values(end),0));
%         sq.find('Bias N/S').after(t_off,sq.minjerk(t_off,sq.find('bias N/S').values(end),0));
%         sq.find('Bias E/W').after(t_off,sq.minjerk(t_off,sq.find('bias E/W').values(end),0));
%     end
% 
%     t_off = linspace(0,40e-3,40);
%     sq.anchor(time_at_evap_end);
%     sq.find('state prep ttl').before(InTrapDelay1,1).after(MW1Duration,0);
%     if opt.mw.analyze(1) ~= 1
%         % % % Turn on repump for in-trap blow away of remaining F = 1 atoms
%         sq.find('Repump Amp TTL').before(InTrapDelay1 - MW1Duration,1).after(BlowDuration1,0);
%         sq.find('Top Repump Shutter').before(InTrapDelay1 + 3e-3 - MW1Duration,0).after(BlowDuration1 + 1e-3,1);
%         sq.find('repump freq').before(InTrapDelay1 - MW1Duration,convert.repump_freq(0));
%         sq.find('Repump Amp').before(InTrapDelay1 - MW1Duration,10);
%     end
% end
% if opt.mw.enable(2) == 1 && opt.mw.analyze(1) ~= 1 % % % Transfer |2,0> -> |1,0>
%     sq.anchor(time_at_evap_end);
%     sq.find('R&S list step trig').before(InTrapDelay2+40e-3,0);
%     sq.find('state prep ttl').before(InTrapDelay2,1).after(MW2Duration,0);
%     if opt.mw.analyze(2) ~= 1
%         % % % Turn on trap for in-trap blow away of remaining F = 2 atoms
%         sq.find('3D MOT Amp TTL').before(InTrapDelay2 - MW2Duration,1).after(BlowDuration2,0);
%         sq.find('3D MOT Amp').before(InTrapDelay2 - MW2Duration,5).after(BlowDuration2,0);
%         sq.find('3D MOT Freq').before(InTrapDelay2 - MW2Duration,RunConversions.mot_freq(0)).after(BlowDuration2,RunConversions.mot_freq(-75));
%     end
% end

%% out of trap MW transfer
% % % Inputs
% pulse 1
MW1Duration = 700*1e-6;
Delay1 = 5e-3;
BlowDuration1 = 1.0*1e-3;

BiasEW = 4;
BiasUD = RunConversions.UD_A_to_V(13);
BiasNS = 0;
MWBiasRampTime = 50e-3;
MWBiasDelay = 100e-3; %100


% pulse 2
MW2Duration = 210*1e-6;
Delay2 = 40e-3;
BlowDuration2 = 100*1e-6;

if opt.mw.enable(1) == 1 % % % Transfer |1,-1> -> |2,0>    % Bias
    sq.anchor(time_at_evap_end + Delay1 - MWBiasDelay);
    t_on = linspace(0,MWBiasRampTime,100);
    t_off = linspace(0,30e-3,20);
    sq.find('Bias U/D').after(t_on,sq.minjerk(t_on,sq.find('bias U/D').values(end),BiasUD));
    sq.find('Bias N/S').after(t_on,sq.minjerk(t_on,sq.find('bias N/S').values(end),BiasNS));
    sq.find('Bias E/W').after(t_on,sq.minjerk(t_on,sq.find('bias E/W').values(end),BiasEW));
    sq.anchor(time_at_evap_end + Delay1 + MW1Duration + 100e-6);
    sq.find('Bias U/D').after(t_off,sq.minjerk(t_off,sq.find('bias U/D').values(end),RunConversions.UD_A_to_V(13)));
    sq.find('Bias N/S').after(t_off,sq.minjerk(t_off,sq.find('bias N/S').values(end),0));
    sq.find('Bias E/W').after(t_off,sq.minjerk(t_off,sq.find('bias E/W').values(end),0));


    t_off = linspace(0,40e-3,40);
    sq.anchor(time_at_evap_end + Delay1);
    sq.find('state prep ttl').set(1);
    sq.delay(MW1Duration);
    sq.find('state prep ttl').set(0);
    sq.find('R&S list step trig').set(0);
    
    if opt.mw.analyze(1) ~= 1
        % % % Turn on repump for in-trap blow away of remaining F = 1 atoms
        ShutterDelay = 3e-3;
        sq.delay(-ShutterDelay);
        sq.find('Top Repump Shutter').set(0);
        sq.find('repump freq').set(convert.repump_freq(0));
        sq.find('Repump Amp').set(10);
        sq.delay(ShutterDelay);
        sq.find('Repump Amp TTL').set(1);
        sq.delay(BlowDuration1);
        sq.find('Repump Amp TTL').set(0);
        sq.find('Top Repump Shutter').set(1);
        sq.find('Repump Amp').set(0);

%         sq.find('Repump Amp TTL').after(MW1Duration,1).after(BlowDuration1,0);
%         sq.find('Top Repump Shutter').before(ShutterDelay,0).after(ShutterDelay + BlowDuration1 + 1e-3,1);
%         sq.find('repump freq').before(ShutterDelay,convert.repump_freq(0));
%         sq.find('Repump Amp').before(ShutterDelay,10);
    end
end
if opt.mw.enable(2) == 1 && opt.mw.analyze(1) ~= 1 % % % Transfer |2,0> -> |1,0>
    sq.anchor(time_at_evap_end + Delay1 + Delay2);
    sq.find('state prep ttl').set(1).after(MW2Duration,0);
    if opt.mw.analyze(2) ~= 1
        % % % Turn on trap for in-trap blow away of remaining F = 2 atoms
        sq.find('3D MOT Amp TTL').after(MW2Duration,1).after(BlowDuration2,0);
        sq.find('3D MOT Amp').set(5).after(MW2Duration + BlowDuration2,0);
        sq.find('3D MOT Freq').set(RunConversions.mot_freq(0)).after(MW2Duration + BlowDuration2,RunConversions.mot_freq(-75));
    end
end

% Stern-Gerlach pulse to test MW transfer
if opt.mw.enable_sg == 1
    if strcmpi('In-Trap',opt.misc.DropCamera) == 1
        % If imaging using the "in-trap" location, use a short SG delay
        SGDelay = 5e-3;
        SGDuration = 0.5e-3;
        Amp = 0.1;

        sq.anchor(timeAtDrop + SGDelay);
        sq.find('mot coil ttl').set(1);
        t = linspace(0,SGDuration,40);
        sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,Amp)));
        sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
        sq.delay(2*SGDuration);
        sq.find('mot coil ttl').set(0);
        sq.find('3d coils').set(convert.mot_coil(0));
    else
        % If not using the in-trap location, use a long SG delay
        SGDelay = 50e-3;
        SGDuration = 0.75e-3;
        mot_off = 0;
        sq.anchor(timeAtDrop + SGDelay);
        sq.find('mot coil ttl').set(1);
        t = linspace(0,SGDuration,40);
        sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,mot_off,mot_off)));%5
        sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(mot_off)));
        sq.delay(2*SGDuration);
        sq.find('mot coil ttl').set(0);
        sq.find('3d coils').set(convert.mot_coil(mot_off));
    end
end


%% Bragg shit flag
BraggOnOff = opt.raman;

k = 2*pi*384229441689483/const.c;  %Frequency of Rb-85 F=3 -> F'=4 transition
Ch1Max = 233;
Ch2Max = 232;

P2onP1 = 1;
P_total = 50;
P1 = P_total/(1+P2onP1);
P2 = (P_total*P2onP1)/(1+P2onP1);
AOMSetting = P_total/Ch2Max;

tau = 30e-6;
T = 1e-3;
T_Sep = 1e-3;
t0 = 40e-3;

braggOrder = 1;
chirp = 2.5106258428e7;

if BraggOnOff == 1
    sq.anchor(timeAtDrop);
    if opt.StatePrep == 0
        sq.find('DDS Trig').set(0).after(10e-3,1);
        sq.ddsTrigDelay = timeAtDrop;
    end

    makeBraggSequence(sq.dds,'k',k,'dt',1e-6,'t0',t0,'T',T,...
        'width',opt.bragg.width,'Tasym',0,'phase',[0,0,opt.bragg.phase],'chirp',chirp,...
        'power',AOMSetting*[1,0,0],'order',braggOrder);
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
    if strcmpi('In-Trap',opt.misc.DropCamera) == 0
        makeImagingSequence(sq,'type','drop 2','tof',opt.tof,...
            'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
            'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.2,...
            'manifold',1,'includeDarkImage',true,'cycle time',150e-3);
    else
        makeImagingSequence(sq,'type','in-trap','tof',opt.tof,...
            'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
            'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.2,...
            'manifold',2,'includeDarkImage',true,'cycle time',150e-3);
    end
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