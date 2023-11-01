function sq = makeBEC_rhys_backup(varargin)


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

%% Stage selection
MOT_status = opt.MOT_status;
CMOT_status = opt.CMOT_status;
PGC_status = opt.PGC_status;
LoadMagTrap_status = opt.LoadMagTrap_status;
MagEvaporation_status = opt.MagEvaporation_status;
LoadOpticalTrap_status = opt.LoadOpticalTrap_status;
OpticalEvaporation_status = opt.OpticalEvaporation_status;
BECCompression_status = opt.BECCompression_status;
MagneticInsensitive_status = opt.MagneticInsensitive_status;

%% Create a conversion object to handle conversions to volts
convert = RunConversions;
imageVoltage = convert.imaging(opt.detuning - 1.8);

%% Initialize sequence - defaults should be handled here
sq = initSequence;


UD = @(x) x* -0.2031 + 2.707; %this converts the input value in amps to the appropriate voltage
interferometry = 0;
DipoleLoad = 1;
RamanF1PumpTest = 0;
RamanImagingSequence = 0;


%% MOT
% % % Turn on dipoles at start of run
sq.find('Imaging Freq').set(convert.imaging(opt.detuning));
sq.find('50w ttl').set(1);
sq.find('25w ttl').set(1);
sq.find('50w amp').set(convert.dipole50(22));
sq.find('25w amp').set(convert.dipole25(17));
sq.find('Top Repump Shutter').set(0);


% % % Set 2D MOT Parameters
% Digital
sq.find('2D MOT Amp TTL').set(1);
sq.find('Push Amp TTL').set(1);
% Analogue
sq.find('Push Freq').set(9.5); % 9.5
sq.find('2D MOT Freq').set(7.65); %7.65

% % Set 3D Mot Parameters
% Digital
sq.find('3D MOT Amp TTL').set(1);
sq.find('Repump Amp TTL').set(1);
sq.find('MOT Coil TTL').set(1);
% Analogue
sq.find('3d coils').set(convert.mot_coil(1.9));
sq.find('3D MOT Amp').set(5);
sq.find('3D MOT Freq').set(convert.mot_freq(-17));
sq.find('Repump Amp').set(9); %9
sq.find('Repump Freq').set(convert.repump_freq(-1.25)); %1.5
sq.find('Liquid Crystal Repump').set(7);

% % % Mag Biases
sq.find('Bias E/W').set(0);
sq.find('Bias N/S').set(0);
sq.find('Bias U/D').set(0);


% % % Duration
sq.delay(opt.MOT_LoadTime);


%% Compressed MOT stage
if CMOT_status == 1
    % % % Goal: produce a smaller/denser cloud

    % stop blowing hot atoms into the MOT
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);
    % Reduce re-radiation pressure (i.e. detune beams)
    sq.find('3D MOT freq').set(convert.mot_freq(-23)); %-25
    sq.find('repump freq').set(convert.repump_freq(-7));
    % Increase/Decrease mag field to compress trapped atoms
    sq.find('3D coils').set(convert.mot_coil(1.45));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(2);
    sq.find('bias u/d').set(UD(0)); %13.4285
    % Duration
    Tcmot = 10.5e-3;
    sq.delay(Tcmot);
end

%% PGC stage
if PGC_status == 1
    % Duration/Ramp
    Tpgc = 37*1e-3;
    t = linspace(0,Tpgc,50);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    % Detune/reduce light
    sq.find('3D MOT Amp').after(t,f(5,3.8));
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-72)));
    sq.find('Repump Amp').set(7);
    sq.find('repump freq').set(convert.repump_freq(-9));

    % Have B = 0 at the atoms
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),0));
    sq.find('bias u/d').set(UD(0));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);

    sq.delay(Tpgc);
end


if RamanF1PumpTest == 1
    % Pump atoms into the F = 1 state
    % i.e. turn repump off
    T = 5*1e-3;
    sq.find('3d coils').set(convert.mot_coil(1.9));
    sq.find('repump amp ttl').set(0);
    sq.find('liquid crystal repump').set(-2.22); %2.22
    sq.find('3D MOT Freq').set(convert.mot_freq(-18));
    sq.find('3D MOT Amp').set(5);
    sq.delay(T);

    %     % Pump atoms into the F = 2 state
    %          % i.e. turn trapping off
    %     T = 30e-3;
    %     sq.find('3D MOT Amp TTL').set(0);
    %     sq.find('repump freq').set(convert.repump_freq(0));
    %     sq.find('Repump Amp').set(9);
    %     sq.delay(T);

end

%% Load into magnetic trap
if LoadMagTrap_status == 1
    % % % De-pump
    %     T = 1e-3;
    %     sq.find('repump amp ttl').set(0);
    %     sq.find('liquid crystal repump').set(-2.22); %2.22
    %     sq.delay(T);

    Tdepump = 1e-3;
    sq.find('repump amp ttl').set(0);
    sq.find('Top repump shutter').set(1);
    sq.find('liquid crystal repump').set(-2.22);

    sq.find('bias u/d').set(UD(0));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);
    sq.delay(Tdepump);


    % % % Mag Load
    sq.find('liquid crystal bragg').set(-3.9);
    sq.find('3D mot amp ttl').set(0);
    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').set(2.05); %2.05

    if MagEvaporation_status == 0
        sq.delay(100e-3);
    end
end

%% Microwave evaporation
if MagEvaporation_status ==1
    evapRate = 10;
    InitialFreq =32; %32
    FinalFreq = 5; % 5
    Tevap = (InitialFreq-FinalFreq)/evapRate;

    t = linspace(0,Tevap,100);
    sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,InitialFreq,FinalFreq)));
    sq.find('mw amp ttl').set(1);
    sq.delay(Tevap);
end

%% Weaken trap while MW frequency fixed
if LoadOpticalTrap_status == 1
    Trampcoils = 180e-3;
    t = linspace(0,Trampcoils,50);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),convert.mot_coil(1)));
    sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
    sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),-0.12*0));
    sq.delay(Trampcoils);
end

%% Optical evaporation
if OpticalEvaporation_status == 1
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
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(4.8),0.5));
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(4.8),0.5));

    sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),10));
    sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
    sq.find('bias u/d').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias u/d').values(end),0));
    sq.delay(Tevap);

    % sq.delay(opt.params(1)); %test the lifetime in the optical trap
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

%% Interfeormetry
if interferometry == 1

    sq.anchor(timeAtDrop);

    sq.find('dds trig').before(10e-3,1);
    sq.find('dds trig').after(10e-3,0); %MOGLabs DDS triggers on falling edge
    sq.find('dds trig').after(10e-3,1);
    sq.ddsTrigDelay = timeAtDrop;

    k = 2*pi*384.2292416894837e12/const.c;  %Frequency of Rb-85 F=3 -> F'=4 transition

    Drop_Time = opt.tof;

    vrel = 2*const.hbar*k/const.mRb;
    chirp = 25.106258428e6;
    T_int = 1e-3;
    Tsep = 27.5e-3; % T sep > sig/v = (sig*m)/(hbar*k), sig = cloud size

    tvs = Drop_Time - 2*T_int - Tsep;   %Velocity pulse time before the Bragg pulse
    %     tvs = 0;   %Velocity pulse time before the Bragg pulse

    Tasym =0;



    %     makeVelocitySelectionpulse(sq.dds,'k',k,'tvs',t0+0e-3,'width',60e-6,'order',2,'chirp',chirp,'power',0.46);

    % % % Add Intensity Noise
    %     I_factor = MakeIntensityNoise('type','acceleration','amp',10,'width',30e-6,'dt',1e-6,'NumPulseWidths',2);
    % % % Add Intensity Ramp
    I_factor = MakeIntensityRamp('amp',10,'width',30e-6,'dt',1e-6,'NumPulseWidths',2);
    % I_factor = 1;

    MakePulseSequence_Rhys(sq.dds,'k',k,'t0',tvs,'T',T_int,'width',30e-6,...
        'Tasym',Tasym,'phase',[0,0,0],'start_order',0,'order',1,'chirp',chirp,...
        'power',0.2*[1,1,0],'PulseType','Square','NumPulseWidths',2,'I_factor',I_factor,'I_ratio',1.5);

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
            'tof1',opt.tof,'repump delay',1e-3,'tof2',8e-3,'cycle time',150e-3,...
            'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',50e-6,...
            'imaging freq',imageVoltage,'repump freq',4.3,'includeDarkImage',true);
    else
        makeImagingSequence(sq,'type',Abs_Analysis_parameters.camera,'tof',opt.tof,...
            'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
            'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.3,...
            'manifold',1,'includeDarkImage',true,'cycle time',150e-3);
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