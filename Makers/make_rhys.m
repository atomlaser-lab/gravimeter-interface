function sq = makeBEC_rhys(varargin)

%% initialise
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
sq.find('Imaging Freq').set(convert.imaging(opt.detuning));


%% Set up the MOT loading values
% Turn Dipoles On At Start of Run
sq.find('50w ttl').set(1);
sq.find('25w ttl').set(1);
sq.find('50w amp').set(convert.dipole50(22));
sq.find('25w amp').set(convert.dipole25(17));



% % Set 2D MOT Parameters
% Digital
sq.find('2D MOT Amp TTL').set(1);
sq.find('Push Amp TTL').set(1);
% Analogue
sq.find('Push Freq').set(9.6);
sq.find('2D MOT Freq').set(7.9);

% % Set 3D Mot Parameters
% Digital
sq.find('3D MOT Amp TTL').set(1);
sq.find('Repump Amp TTL').set(1);
sq.find('MOT Coil TTL').set(1);
% Analogue
sq.find('3d coils').set(convert.mot_coil(1.7));
sq.find('3D MOT Amp').set(5); 
sq.find('3D MOT Freq').set(convert.mot_freq(-19));
sq.find('Repump Amp').set(9);
sq.find('Repump Freq').set(convert.repump_freq(-1));
sq.find('Liquid Crystal Repump').set(7);

% % % Biases
sq.find('Bias E/W').set(0*3); 
sq.find('Bias N/S').set(0); %optimum in opposite direction
sq.find('Bias U/D').set(0); %optimum in opposite direction

% % % Duration
sq.delay(opt.MOT_LoadTime);

%% Compressed MOT stage
if CMOT_status == 1
    % % % Goal: produce a smaller/denser cloud
    % stop blowing hot atoms into the MOT
    sq.find('2D MOT Amp TTL').before(5e-3,0);
    sq.find('Push Amp TTL').before(5e-3,0);

    % Reduce re-radiation pressure (i.e. detune beams) 
    sq.find('3D MOT freq').set(convert.mot_freq(-40)); %-25?
    sq.find('repump freq').set(convert.repump_freq(-6)); 

    % Increase mag field to compress trapped atoms
    sq.find('3d coils').set(convert.mot_coil(1.5));

    % Duration
    TCMOT = 12e-3;
    sq.delay(TCMOT);
end

%% PGC stage
if PGC_status == 1
    % % % Goal: Produce a low temperature cloud
    % Duration
    Tpgc = 30e-3;
    % Define a ramp
    t = linspace(0,Tpgc,100);
    f = @(vi,vf) sq.linramp(t,vi,vf);    

    % Turn mag fields off
    sq.find('Bias E/W').set(0);
    sq.find('Bias N/S').set(0);
    sq.find('Bias U/D').set(0);

    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),convert.mot_coil(0.525))); % 0.5 - 0.525 - 0.55
    
    % Detune light
    sq.find('3D MOT Amp').after(t,f(5,3)); %3.5 is optimum for PGC, 3 is optimum for mag evap
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-65))); %-60 for PGC, -65 for mag evap
    sq.find('repump freq').set(convert.repump_freq(-8.5)); %-9, -8.5 for mag evap   
    sq.delay(Tpgc);
end

%% Load into magnetic trap
if LoadMagTrap_status == 1
    % De-pump the atoms (put into F = 1, mf = -1) by turning off repump
    T_depump = 1e-3;
    sq.find('Repump Amp TTL').set(0);
    sq.find('Liquid Crystal Repump').set(-2.2);
    sq.delay(T_depump);

    % Turn off trapping light while in mag trap
    sq.find('3D MOT Amp TTL').set(0);
    sq.find('Repump Amp TTL').set(0);
    sq.find('liquid crystal bragg').set(-3.9);
    

    % Define mag ramp
    t = linspace(0,50e-3,100);
    f = @(vi,vf) sq.linramp(t,vi,vf); 
    StartCurrent = 7;
    EndCurrent = opt.params; %13 want low field for mag load, but high field for rethermalisation
    sq.find('3D coils').after(t,f(convert.mot_coil(StartCurrent),convert.mot_coil(EndCurrent)));
    if MagEvaporation_status == 0
        sq.delay(120e-3);
    end
end

%% Microwave evaporation
if MagEvaporation_status ==1

end

%% Weaken trap while MW frequency fixed
if LoadOpticalTrap_status == 1

end

%% Optical evaporation
if OpticalEvaporation_status == 1

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
interferometry = 0;
if interferometry == 1
    k = 2*pi*384.2292416894837e12/const.c;  %Frequency of Rb-85 F=3 -> F'=4 transition

    Drop_Time = 10e-3;
    TOF = opt.tof;

    vrel = 2*const.hbar*k/const.mRb;
    dv = 700e-6/TOF;
    dv = 700e-6/216.5e-3;
    chirp = 25.106258428e6;
    T_int =10e-3;
    tvs = Drop_Time - 2*T - Tsep;   %Velocity pulse time before the Bragg pulse
    Tasym =0;

    sq.dds.anchor(timeAtDrop);

    %     makeVelocitySelectionpulse(sq.dds,'k',k,'tvs',t0+0e-3,'width',60e-6,'order',2,'chirp',chirp,'power',0.46);

    makeBraggSequence(sq.dds,'k',k,'t0',tvs,'T',T_int,'width',40e-6,...
        'Tasym',Tasym,'phase',[0,0,0],'start_order',0,'order',1,'chirp',chirp,...
        'power',0.2*[1,0,0]);
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
    makeImagingSequence(sq,'type',Abs_Analysis_parameters.camera,'tof',opt.tof,...
        'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
        'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.3,...
        'manifold',1,'includeDarkImage',true,'cycle time',150e-3);
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