function sq = makeBEC_Bragg_NewCloud(varargin)
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

sq.find('3D MOT Freq').set(convert.mot_freq(-18.5 + 0.5));
sq.find('3D MOT Amp').set(RunConversions.mot_power(1 - 0.4));
sq.find('Repump freq').set(convert.repump_freq(-1.75 + 0.25));
sq.find('Repump Amp').set(RunConversions.repump_power(1 - 0.1));

sq.find('MOT coil TTL').set(1);
sq.find('3d coils').set(RunConversions.mot_coil(1.35));
sq.find('bias u/d').set(convert.UD_A_to_V(13.8 - 0.5 -0.25));
sq.find('bias e/w').set(2 + 0.5 -1.25);
sq.find('bias n/s').set(0);

sq.delay(opt.MOT_LoadTime);

%% Compressed MOT stage
if opt.CMOT_status == 1
    PushDelay = 2*1e-3;
    sq.find('2D MOT Amp TTL').before(PushDelay,0);
    sq.find('push amp ttl').before(PushDelay,0);

    sq.find('3D MOT freq').set(convert.mot_freq(-23.5));
    sq.find('3D MOT Amp').set(RunConversions.mot_power(0.8 - 0.1));
    sq.find('repump freq').set(convert.repump_freq(-9 + 0.5));
    sq.find('Repump Amp').set(RunConversions.repump_power(1 - 0.1));
    
    sq.find('3D coils').set(RunConversions.mot_coil(1.5));
    sq.find('bias e/w').set(1 + 2);
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

    sq.find('3D MOT Amp').after(t,f(RunConversions.mot_power(1),RunConversions.mot_power(0.8))); %+ 0.2
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-65 - 1)));
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),RunConversions.mot_coil(0.6 - 0.05)));
    sq.find('repump freq').set(convert.repump_freq(-9 - 0.25 + 0.2));

    sq.find('bias u/d').set(convert.UD_A_to_V(11));
    sq.find('bias e/w').set(1 + 2);
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
    evapRate = 12 - 1;
    evapStart = 34;
    evapEnd = 8 - 3;
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
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(opt.dipoles),0.56)); %1.47
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(opt.dipoles),0.72)); 
    sq.delay(Tevap);
    time_at_evap_end = sq.time;    
end

%% Drop atoms
timeAtDrop = sq.time;
sq.anchor(timeAtDrop);

sq.find('3D mot amp ttl').set(0);
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
%% out of trap MW transfer
% % % Inputs
% pulse 1
MW1Duration = 500e-6;
% MW1Duration = 700*1e-6;
Delay1 = 11*1e-3;

BlowDuration1 = 0.05*1e-3; % 1

BiasEW = 4;
BiasUD = RunConversions.UD_A_to_V(13);
BiasNS = 0;
MWBiasRampTime = 50e-3;
MWBiasDelay = 100e-3; %100


% pulse 2
MW2Duration = 210*1e-6;
Delay2 = (32 + 2.725)*1e-3; %32
BlowDuration2 = 5*1e-6; %100

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
    end
end
if opt.mw.enable(2) == 1 && opt.mw.analyze(1) ~= 1 % % % Transfer |2,0> -> |1,0>
    sq.anchor(time_at_evap_end + Delay2);
    sq.find('state prep ttl').set(1).after(MW2Duration,0);
    if opt.mw.analyze(2) ~= 1
        % % % Turn on trap for in-trap blow away of remaining F = 2 atoms
        sq.find('3D MOT Amp TTL').after(MW2Duration,1).after(BlowDuration2,0);
        sq.find('3D MOT Amp').set(5).after(MW2Duration + BlowDuration2,0);
        sq.find('3D MOT Freq').set(RunConversions.mot_freq(-5)).after(MW2Duration + BlowDuration2,RunConversions.mot_freq(-75));
    end
end

if Delay2 - Delay1 < 21e-3
    warning('Insufficient delay between MW pulse 1 and 2. MKU doesnt have enough time to change isntructions')
end
% Stern-Gerlach pulse to test MW transfer
if opt.mw.enable_sg == 1
    if strcmpi('In-Trap',opt.misc.DropCamera) == 1
        SGDelay = 10e-3;
        SGDuration = 4e-3;
        Amp = 1; %0.5
        
        sq.anchor(timeAtDrop + SGDelay);
        sq.find('mot coil ttl').set(1);
        t = linspace(0,SGDuration,40);
        sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,Amp)));%5
        sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
        sq.delay(2*SGDuration);
        sq.find('mot coil ttl').set(0);
        sq.find('3d coils').set(convert.mot_coil(0));

    else
        % If not using the in-trap location, use a long SG delay
        SGDelay = 50e-3;
        SGDuration = 0.75e-3;
        mot_on = 0.05;        
        mot_off = 0;
        sq.anchor(timeAtDrop + SGDelay);
        sq.find('mot coil ttl').set(1);
        t = linspace(0,SGDuration,40);
        sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,mot_on)));%5
        sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(mot_off)));
        sq.delay(2*SGDuration);
        sq.find('mot coil ttl').set(0);
        sq.find('3d coils').set(convert.mot_coil(mot_off));
    end

    if opt.mw.enable(2) == 0 && opt.mw.enable(1) == 1
        if SGDelay < Delay1
            warning('SG pulse occurs during MW pulse');
        end
    elseif opt.mw.enable(2) == 1 && opt.mw.enable(1) == 1
        if SGDelay < Delay1 + Delay2
            warning('SG pulse occurs during MW pulse');
        end
    end
end


%% Bragg flag
k = 2*pi*384229441689483/const.c;  %Frequency of Rb-85 F=3 -> F'=4 transition
Ch1Max = 215;
Ch2Max = 240;

% P_total = 112*opt.params(1);
% FreqError = opt.params(2); %in Hz, max ~ 3*(4wr)/4 ~ 3 wr
% P_total = 112.5;
% P_total = opt.params;

% P_total = opt.params;
FreqError = 0;
% FreqError = opt.params; %in Hz, max ~ 3*(4wr)/4 ~ 3 wr
% AOMSetting = P_total/Ch2Max;

P2onP1 = 1;
P1 = P_total/(1+P2onP1);
P2 = (P_total*P2onP1)/(1+P2onP1);
AOM1 = P1/Ch2Max;
AOM2 = P2/Ch2Max;

if P1 > Ch1Max || P2 > Ch2Max
    warning('Not enough power')
    return
end

tau = 30*1e-6;
T = 10e-3;
t0 = 35e-3; %35
braggOrder = 1;
PulseType = PulseTypes.Primitive;
NumStd = 1.5;

chirp = 2.5106258428e7; %used in RT's makesequence
% chirp = 2.511010143833653*1e7; %found on 27/08/2024
chirp = 2.511010323051248*1e7; %found on 24/09/2024
% chirp = opt.params*1e7;

if opt.raman == 1
    sq.anchor(timeAtDrop);
    if opt.StatePrep == 0
        sq.find('DDS Trig').set(0).after(10e-3,1);
        sq.ddsTrigDelay = timeAtDrop;
    end

%     makeBraggCompositePulse_Rhys6(sq.dds,'k',k,'chirp',chirp,...
%         'freqerror',FreqError,'pulsetype',PulseType,'numstd',NumStd,'pulseShape','Gaussian', ...
%         'dt',1e-6,'t0',t0,'T',T,'tasym',0,'width',tau,...
%         'phase',0,'power1',AOM1*1,'power2',AOM2*1,'order',braggOrder,...
%         'NoiseType','acceleration','normdisp',0,'BeamRadius',10e-3,'t0_effective',7.5e-3,'RampOnOff',0);

%     % accelamp or normdisp
%     makeBraggSequence_Rhys2(sq.dds,'k',k,'chirp',chirp,...
%         'freqerror',FreqError,...
%         'dt',1e-6,'t0',t0,'T',T,'tasym',0,'width',tau,...
%         'phase',[0,0,0],'power1',AOM1*[1,0,0],'power2',AOM2*[1,0,0],'order',braggOrder,...
%         'NoiseType','acceleration','normdisp',0,'BeamRadius',10e-3,'t0_effective',7.5e-3,'RampOnOff',0);
% 
    makeBraggSequence_Rhys2(sq.dds,'k',k,'chirp',chirp,...
        'dt',1e-6,'t0',t0,'T',T,'tasym',0.e-3,'width',tau,'freqerror',FreqError,...
        'phase',[0,0,opt.params],'power1',AOM1*[0.5,1,0.5],'power2',AOM2*[0.5,1,0.5],'order',braggOrder,...
        'NoiseType','acceleration','normdisp',0,'BeamRadius',10e-3,'t0_effective',7.5e-3,'RampOnOff',0);
% % %     makeBraggSequence_Rhys2(sq.dds,'k',k,'chirp',chirp,...
% % %         'dt',1e-6,'t0',t0,'T',T,'tasym',0.e-3,'width',tau,'freqerror',FreqError,...
% % %         'phase',[0,0,0],'power1',AOM1*[1,0,0],'power2',AOM2*[1,0,0],'order',braggOrder,...
% % %         'NoiseType','acceleration','normdisp',0,'BeamRadius',10e-3,'t0_effective',7.5e-3,'RampOnOff',0);

%     makeBraggSequence_Rhys2(sq.dds,'k',k,'chirp',chirp,...
%         'dt',1e-6,'t0',t0,'T',T,'tasym',0,'width',tau,...
%         'phase',[0,0,opt.params(1)],'power1',AOM1*[0.5,1,0.5],'power2',AOM2*[0.5,1,0.5],'order',braggOrder,...
%         'NoiseType','acceleration','normdisp',opt.params(2),'BeamRadius',2.5e-3,'t0_effective',7.5e-3,'RampOnOff',1);
%   
%     makeBraggSequence_Rhys2(sq.dds,'k',k,'chirp',chirp,...
%         'dt',1e-6,'t0',t0,'T',T,'tasym',0,'width',tau,...
%         'phase',[0,0,90],'power1',AOM1*[0.5,1,0.5],'power2',AOM2*[0.5,1,0.5],'order',braggOrder,...
%         'NoiseType','acceleration','normdisp',1.2,'BeamRadius',2.5e-3,'t0_effective',7.5e-3,'RampOnOff',1);

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
            'manifold',1,'includeDarkImage',true,'cycle time',150e-3);
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