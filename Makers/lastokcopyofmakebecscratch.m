function sq = makeBEC_scratch(varargin)
%% Ensure SequenceOptions are available
% opt = ensureSequenceOptionsAvailable();

allVarNames = evalin('base', 'who');

% Check if any variable is an instance of SequenceOptions
seqOptVarIndices = cellfun(@(varName) isa(evalin('base', varName), 'SequenceOptions'), allVarNames);
seqOptVarNames = allVarNames(seqOptVarIndices);

if isempty(seqOptVarNames)
    if nargin == 1 && ~exist(varargin{1},'var')
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

%Test selection
LoadMagTrap_test = 0;
dipoleposition_test = 0;

%% Create a conversion object to handle conversions to volts
convert = RunConversions;
% imageVoltage = convert.imaging(opt.detuning - 1.8);
imageVoltage = convert.imaging(opt.detuning);

UD = @(x) x* -0.2031 + 2.707; %this converts the input value in amps to the appropriate voltage
% UDreverse = @(y) y*-4.924 + 13.33; %This is just the inverted function to get the amps out of the voltage in case needed.
%U/D bias field converter to amps to volts (possible values are -0.58 A to 14 A) which corresponds to a voltage from (2.823V to -0.14V)

%% Initialize sequence - defaults should be handled here
sq = initSequence;

if MOT_status == 1

    sq.find('Top repump shutter').set(0);
    sq.find('2D MOT Freq').set(7.705);
    sq.find('Push Freq').set(9.45);
    sq.find('liquid crystal repump').set(7);
    sq.find('Imaging Freq').set(convert.imaging(opt.detuning));
    sq.find('3D MOT Freq').set(convert.mot_freq(-19));    %Use -25 MHz for 4 s loading times, -27.5 MHz for 6 s
    sq.find('Repump Amp').set(9);
    sq.find('Repump freq').set(convert.repump_freq(-1.625));
    sq.find('50w ttl').set(0);
    sq.find('25w ttl').set(0);
%     sq.find('50w amp').set(convert.dipole50(27));
%     sq.find('25w amp').set(convert.dipole25(23.5));


    %% Set up the MOT loading values

    sq.find('MOT coil TTL').set(1);     %Turn on the MOT coils
    sq.find('3d coils').set(convert.mot_coil(2.2)); %2.2

    %     sq.find('bias u/d').set(UD(13.5)); %13.4285
    sq.find('bias u/d').set(0); %13.4285
    sq.find('bias e/w').set(1);
    sq.find('bias n/s').set(2);

    Tmot = opt.MOT_LoadTime;                           %6 s MOT loading time
    sq.delay(Tmot);                     %Wait for Tmot
end
%% Compressed MOT stage

if CMOT_status == 1
    % %Turn off the 2D MOT and push beam 10 ms before the CMOT stage
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);

    %Increase the cooling and repump detunings to reduce re-radiation
    %pressure, and weaken the trap
    sq.find('3D MOT freq').set(convert.mot_freq(-23)); %-25
    sq.find('repump freq').set(convert.repump_freq(-7));
    sq.find('3D coils').set(convert.mot_coil(1.45));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(2);
%     sq.find('bias u/d').set(UD(0)); %13.4285
    sq.find('bias u/d').set(0); %13.4285

    %
    Tcmot = 10.5e-3;                      %12.5 ms CMOT stage
    sq.delay(Tcmot);                    %Wait for time Tcmot

end

%% PGC stage
if PGC_status == 1

    Tpgc = 37*1e-3;
    %Define a function giving a 100 point smoothly varying curve
    t = linspace(0,Tpgc,50);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    %Smooth ramps for these parameters
    sq.find('3D MOT Amp').after(t,f(5,3.8));
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-72)));
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),0));
    sq.find('repump freq').set(convert.repump_freq(-9.5)); % was 9 and gives a slight i
    sq.find('Repump Amp').set(7);

    sq.find('bias u/d').set(0);
%     sq.find('bias u/d').set(UD(0));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);

    sq.delay(Tpgc);

%     %      Turn off the repump field for optical pumping - 1 ms
    Tdepump = 1e-3;
    sq.find('repump amp ttl').set(0);
%     sq.find('Top repump shutter').set(1);
    sq.find('liquid crystal repump').set(-2.22);
    
    sq.find('bias u/d').set(0);
%     sq.find('bias u/d').set(UD(0));
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);
    sq.delay(Tdepump);

end %end PGC


%% Load into magnetic trap
if LoadMagTrap_status == 1
    sq.find('liquid crystal bragg').set(-3);
    sq.find('3D mot amp ttl').set(0);
    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').set(.7);
    load_ramp = 5e-3;
    t = linspace(0,load_ramp,100);
    f = @(vi,vf) sq.linramp(t,vi,vf);
    sq.find('3D coils').after(t,f(.7,1.4));

    sq.find('MOT coil ttl').set(1);

    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);
%     sq.find('bias u/d').set(UD(0));
    sq.find('bias u/d').set(0);


    sq.delay(5*1e-3);

    if LoadMagTrap_test == 1
        sq.delay(opt.params); %to test lifetime in magtrap
    end

end

%% Microwave evaporation
if MagEvaporation_status ==1 % Provide etunings in MHz from the Rb hyperfine splitting
    evapRate = 16;
    evapStart = 40;
    evapEnd = 5;
    Tevap = (evapStart-evapEnd)/evapRate;
    t = linspace(0,Tevap,120);
    sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,evapStart,evapEnd)));
    sq.find('mw amp ttl').set(1);
    sq.delay(Tevap);
    sq.find('mw amp ttl').set(0);

end %end MagEvaporation

if LoadOpticalTrap_status == 1 % Weaken trap while MW frequency fixed
if dipoleposition_test ==0
    Trampcoils = 180e-3;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),0));
    sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),00));
    sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),00));
    sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),0));
%     sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),UD(14)));

    sq.delay(Trampcoils);

  
elseif dipoleposition_test == 1
    Trampcoils = 180e-3;
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),0.6));
    sq.delay(Trampcoils);
    sq.find('3d coils').set(-.075);
    sq.delay(27.5e-3);
end

end
%% Optical evaporation
if OpticalEvaporation_status == 1
sq.delay(30e-3);
Tevap = 3;
t = linspace(0,Tevap,250);
sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(5.57),.4));
sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(4.65),.4));

% sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),0));
% sq.find('bias n/s').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias n/s').values(end),0));
sq.find('bias u/d').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias u/d').values(end),0));
sq.delay(Tevap);

end

%% Drop atoms
% sq.find('mot coil ttl').set(0);
% sq.find('3D Coils').set(convert.mot_coil(0));
% sq.delay(35e-3 - opt.tof);

timeAtDrop = sq.time; %Store the time when the atoms are dropped for later
sq.anchor(timeAtDrop);
sq.find('3D mot amp ttl').set(0);
% sq.find('bias e/w').before(200e-3,0);
% sq.find('bias n/s').before(200e-3,0);
% sq.find('bias u/d').before(200e-3,0);

sq.find('bias e/w').before(10e-3,0);

sq.find('bias n/s').before(10e-3,0);
sq.find('bias u/d').before(10e-3,0);

sq.find('mw amp ttl').set(0);
sq.find('mot coil ttl').set(0);
sq.find('3D Coils').set(convert.mot_coil(0));
sq.find('25w ttl').set(0);
sq.find('50w ttl').set(0);
sq.find('50w amp').set(convert.dipole50(0));
sq.find('25w amp').set(convert.dipole25(0));
sq.find('Top repump shutter').set(0);
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