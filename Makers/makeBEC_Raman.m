function sq = makeBEC_Raman(varargin)


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
imageVoltage = convert.imaging(opt.detuning);

%% Initialize sequence - defaults should be handled here
sq = initSequence;

%    U/D bias field converter to amps to volts (possible values are -0.58 A to 14 A) which corresponds to a voltage from (2.823V to -0.14V)
UD = @(x) x* -0.2031 + 2.707; %this converts the input value in amps to the appropriate voltage
% UDreverse = @(y) y*-4.924 + 13.33; %This is just the inverted function to get the amps out of the voltage in case needed.
sq.find('liquid crystal repump').set(3);
sq.find('Top repump shutter').set(0);

sq.find('Imaging Freq').set(convert.imaging(opt.detuning));
sq.find('3D MOT Freq').set(convert.mot_freq(-16.75));    %Use -17 MHz for 4 s loading times,
sq.find('Repump freq').set(convert.repump_freq(-1.5)); %2 %-1.5

if opt.LoadOpticalTrap_status == 1
    sq.find('50w ttl').set(1);
    sq.find('25w ttl').set(1);
    sq.find('50w amp').set(convert.dipole50(35)); %22 %35
    sq.find('25w amp').set(convert.dipole25(20)); %17 %22 %20
end
%% Set up the MOT loading values

sq.find('MOT coil TTL').set(1);     %Turn on the MOT coils
sq.find('3d coils').set(0.15); %42
sq.find('bias u/d').set(0); %-1
% sq.find('bias e/w').set(4);
sq.find('bias n/s').set(0);

Tmot = opt.MOT_LoadTime;                           %6 s MOT loading time
sq.delay(Tmot);                     %Wait for Tmot

%% Compressed MOT stage

if CMOT_status == 1
    % %Turn off the 2D MOT and push beam 10 ms before the CMOT stage
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);

    %Increase the cooling and repump detunings to reduce re-radiation
    %pressure, and weaken the trap
    sq.find('3D MOT freq').set(convert.mot_freq(-22)); %-18 %-19 %-21
    sq.find('repump freq').set(convert.repump_freq(-8));
    sq.find('3D coils').set(0.17);
    % sq.find('bias e/w').set(4);
    sq.find('bias n/s').set(0); %7
    sq.find('bias u/d').set(0);

    %
    Tcmot = 13*1e-3;                      %12.5 ms CMOT stage
    sq.delay(Tcmot);                    %Wait for time Tcmot

end

%% PGC stage
if PGC_status == 1

    Tpgc = 20*1e-3; %20 ms
    %Define a function giving a 100 point smoothly varying curve
    t = linspace(0,Tpgc,50);
    f = @(vi,vf) sq.linramp(t,vi,vf);

    %Smooth ramps for these parameters
    sq.find('3D MOT Amp').after(t,f(5,3.8));
    sq.find('3D MOT Freq').after(t,f(sq.find('3D MOT Freq').values(end),convert.mot_freq(-69))); %66 for 12s
    sq.find('3D coils').after(t,f(sq.find('3D coils').values(end),0));
    sq.find('repump freq').set(convert.repump_freq(-9));


    sq.find('bias u/d').set(3);
    sq.find('bias e/w').set(0);
    sq.find('bias n/s').set(0);

    sq.delay(Tpgc);
    % Turn off the repump field for optical pumping - 1 ms
    Tdepump = 1e-3;
    sq.find('repump amp ttl').set(0);
    sq.find('Top repump shutter').set(1);
    sq.find('liquid crystal repump').set(-2.22);
    % sq.find('bias u/d').set(3);
    % sq.find('bias e/w').set(0);
    % sq.find('bias n/s').set(0);
    sq.delay(Tdepump);
end %end PGC


%% Load into magnetic trap
if LoadMagTrap_status == 1
    sq.find('liquid crystal bragg').set(-3);
    sq.find('3D MOT Amp').set(0);
    sq.find('3D mot amp ttl').set(0);
    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').set(1.6);
    % sq.delay(3.0); %test lifetime in mag trap
end

%% Microwave evaporation
if MagEvaporation_status ==1
    % Provide detunings in MHz from the Rb hyperfine splitting
    sq.delay(22*1e-3);
    evapRate = 14; %10 %12
    evapStart = 40; %37
    evapEnd = 8;
    Tevap = (evapStart-evapEnd)/evapRate;
    t = linspace(0,Tevap,100);
    sq.find('mw freq').after(t,convert.microwave(sq.linramp(t,evapStart,evapEnd)));
    sq.find('mw amp ttl').set(1);
    sq.delay(Tevap);
    % sq.find('mw amp ttl').set(0);
    % sq.delay(3);

end %end MagEvaporation

%% Weaken trap while MW frequency fixed
if LoadOpticalTrap_status == 1
    dipoleposition_test = 0;

    if dipoleposition_test == 0
        Trampcoils = 180e-3;
        t = linspace(0,Trampcoils,100);
        sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),0.5)); %0.4
        sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
        sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
        sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),0));
        sq.delay(Trampcoils);


    elseif dipoleposition_test == 1
        Trampcoils = 180e-3;
        t = linspace(0,Trampcoils,100);
        sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),0.45));%.45
        sq.delay(Trampcoils);
        sq.find('3d coils').set(-.1);
        %         sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
        %         sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
        %         sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),3));
        sq.delay(22.5e-3);
    end
end




%% Optical evaporation
if OpticalEvaporation_status == 1
    % Ramp down magnetic trap in 1 s

    Trampcoils = 0.9; %.8
    t = linspace(0,Trampcoils,100);
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),0.1)); %0.14
    sq.find('mw amp ttl').anchor(sq.find('3d coils').last).before(100e-3,0);
    sq.find('mot coil ttl').at(sq.find('3d coils').last,0);
    % %
    % % At the same time, start optical evaporation
    % %
    sq.delay(5*1e-3); %10
    Tevap =2.5;
    t = linspace(0,Tevap,100);
    sq.find('50W amp').after(t,sq.expramp(t,sq.find('50w amp').values(end),convert.dipole50(1.5),.42)); %1.6 %1.53
    sq.find('25W amp').after(t,sq.expramp(t,sq.find('25w amp').values(end),convert.dipole25(1.40),.7)); %1.26 %1.28 %1.46

    sq.find('bias e/w').after(t(1:end/2),@(x) sq.linramp(x,sq.find('bias e/w').values(end),0));
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
% sq.find('DDS trig').set(0).after(10e-3,1);
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

sq.find('dds trig').before(10e-3,1);
sq.find('dds trig').after(1e-3,0); %MOGLabs DDS triggers on falling edge
sq.find('dds trig').after(10e-3,1);
sq.ddsTrigDelay = timeAtDrop-9e-3;
%% MW state prep
% sq.find('dds trig').before(10e-3,1);
% sq.find('dds trig').after(1e-3,0); %MOGLabs DDS triggers on falling edge
% sq.find('dds trig').after(10e-3,1);
% sq.ddsTrigDelay = timeAtDrop-9e-3;
if any(opt.mw.enable)
    %
    % Apply a pair of microwave pulses to effect the transfers
    % |F=1,m=-1> -> |F=2,m=0> -> |F=1,m=0>.  The first pulse is applied
    % 10 ms after the atoms are dropped to minimize any possible
    % state-changing collisions.  The "R&S list step trig" skips to the
    % next frequency on the rising edge and resets the list on the
    % falling edge
    %
    if opt.mw.enable(1)
        sq.anchor(timeAtDrop);
        sq.find('bias e/w').set(5);
        sq.find('R&S list step trig').set(1);
        sq.delay(4e-3);
        sq.find('state prep ttl').set(1);
        sq.delay(450e-6);
        sq.find('state prep ttl').set(0);



        sq.find('Repump Amp TTL').set(1).after(1e-3,0);
        sq.find('Top repump shutter').before(10e-3,0);
        sq.find('Liquid Crystal Repump').set(7).after(1e-3,-2.22);
        sq.find('repump freq').set(4.3);
        sq.find('Top repump shutter').after(15e-3,1);
    end

    if opt.mw.enable(2) && ~opt.mw.analyze(1)
        sq.find('Repump Amp TTL').set(1).after(1e-3,0);
        sq.find('Top repump shutter').before(10e-3,0);
        sq.find('Liquid Crystal Repump').set(7).after(1e-3,-2.22);
        sq.find('repump freq').set(4.3);
        sq.find('Top repump shutter').after(15e-3,1);
        %         sq.find('R&S list step trig').set(0);
        %         sq.delay(20e-3);
        %         sq.find('state prep ttl').set(1);
        %         sq.delay(250e-6);
        %         sq.find('state prep ttl').set(0);
        %
        %         sq.find('R&S list step trig').set(1);
        %         sq.find('bias e/w').set(0);
        % sq.find('3D MOT Amp').set(3.8);
        %
        %         sq.find('3D MOT Amp TTL').set(1).after(100e-6,0);
    end
else
    sq.find('bias e/w').at(timeAtDrop,0);
end

if opt.mw.enable_sg || opt.mw.analyze(1)
    %
    % Apply a Stern-Gerlach pulse to separate states based on magnetic
    % moment.  A ramp is used to ensure that the magnetic states
    % adiabatically follow the magnetic field
    %
    sq.anchor(timeAtDrop + 15e-3);

    %     sq.delay(30e-3);
    %     sq.waitFromLatest(5e-3);
    Tsg = 2e-3;
    sq.find('mot coil ttl').set(1);
    t = linspace(0,Tsg,20);
    sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,5.5)));
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.delay(2*Tsg);
    sq.find('mot coil ttl').set(0);
    sq.find('3d coils').set(convert.mot_coil(0));
end


%% Raman
Raman = 1;
if Raman == 1

    % Pulse Parameters
    PulseWidth = 26*1e-3 - 1e-3;
    dt = 200e-6;
    delta = 0.18 + opt.params;

%     t0 = 0.007;
  t0 = 4*1e-3;
%     t0 = 20e-3;
%     sq.find('bias u/d').before(100e-6,opt.params);  

%     sq.find('bias e/w').before(100e-6,0);  
    sq.find('bias n/s').before(100e-6,0);
    sq.find('bias u/d').before(100e-6,0);
% 
% MakePulseSequence_Rhys(sq.dds,'t0',t0,'T',1e-3,'width',PulseWidth,'dt',dt,...
%         'phase',[0,0,0],'delta',delta,...
%         'power1',1*[1,0,0],'power2',1*[1,0,0],'PulseType','Square');

%     makeRamanPulseSequence(sq.dds,'t0',t0,'width',PulseWidth,'dt',dt,'delta',delta);

cleaner_time = t0+PulseWidth+0.1e-3;

%     sq.find('bias e/w').after(cleaner_time+PulseWidth/2,0);  

    sq.find('3D MOT freq').after(cleaner_time,RunConversions.mot_freq(0));
    sq.find('3D MOT Amp TTL').after(cleaner_time,1);
    sq.find('3D MOT Amp').after(cleaner_time,5);
    sq.find('3D MOT Amp TTL').after(100e-6,0);
    sq.find('3D MOT Amp').after(100e-6,0);

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
        'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.2,...
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