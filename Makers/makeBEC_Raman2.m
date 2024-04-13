function sq = makeBEC_Raman2(varargin)

% This Transfers atoms from |2,0> to |1,0>

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

% Bragg Calibration data used in initSequence.
% For now I will simply load/set my own calibration data here. I'll create
% another object later
sq.dds(1).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
sq.dds(2).power_conversion_method = DDSChannel.POWER_CONVERSION_DBM_INTERP;
% calibData = load('RamanAOMData_formatted.mat');
calibData = load('RamanAOMData_11042024');
sq.dds(1).calibrationData = calibData.data_ch1;
sq.dds(2).calibrationData = calibData.data_ch2;

%    U/D bias field converter to amps to volts (possible values are -0.58 A to 14 A) which corresponds to a voltage from (2.823V to -0.14V)
UD = @(x) x* -0.2031 + 2.707; %this converts the input value in amps to the appropriate voltage
% UDreverse = @(y) y*-4.924 + 13.33; %This is just the inverted function to get the amps out of the voltage in case needed.
sq.find('liquid crystal repump').set(3);
sq.find('Top repump shutter').set(0);

sq.find('Imaging Freq').set(convert.imaging(opt.detuning));
sq.find('3D MOT Freq').set(convert.mot_freq(-16.75));    %Use -17 MHz for 4 s loading times,
sq.find('Repump freq').set(convert.repump_freq(-1.5)); %2 %-1.5

if opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
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

if opt.CMOT_status == 1
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
if opt.PGC_status == 1

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
if opt.LoadMagTrap_status == 1 && opt.JustMOT ~= 1
    sq.find('liquid crystal bragg').set(-3);
    sq.find('3D MOT Amp').set(0);
    sq.find('3D mot amp ttl').set(0);
    sq.find('MOT coil ttl').set(1);
    sq.find('3D coils').set(1.6);
    % sq.delay(3.0); %test lifetime in mag trap
end

%% Microwave evaporation
if opt.MagEvaporation_status == 1 && opt.JustMOT ~= 1
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
if opt.LoadOpticalTrap_status == 1 && opt.JustMOT ~= 1
    dipoleposition_test = 0;

    if dipoleposition_test == 0
        Trampcoils = 180e-3;
        t = linspace(0,Trampcoils,100);
        sq.find('3d coils').after(t,sq.minjerk(t,sq.find('3d coils').values(end),0.5)); %0.4
        sq.find('bias e/w').after(t,sq.minjerk(t,sq.find('bias e/w').values(end),0));
        sq.find('bias n/s').after(t,sq.minjerk(t,sq.find('bias n/s').values(end),0));
        sq.find('bias u/d').after(t,sq.minjerk(t,sq.find('bias u/d').values(end),0));
        sq.delay(Trampcoils);

        %Hold as atoms fall out of trap
        sq.delay(120e-3);


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
if opt.OpticalEvaporation_status == 1 && opt.JustMOT ~= 1
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


%% Raman alignment test
RamanAlignment = 0;
if RamanAlignment == 1 && opt.raman.OnOff ~= 1
    Ch2_Pratio = 1;
    Ch1_Pratio = 1;

    % % % Inputs
    % Timing
    TriggerDuration = 1e-3; 
    triggerDelay = 1e-3; 

    PulseWidth = 10e-3; %start large and make smaller as you align
    dt = 1e-3;
    TOF = 15*1e-3;

    % Pulse Parameters    
    chirp = 25.106258428e6;
    k = 22.731334388721734;
    delta = 8; % 12

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
    sq.anchor(timeAtDrop + TOF);
    sq.find('Raman DDS Trig').before(TriggerDuration,1);
    sq.find('Raman DDS Trig').after(TriggerDuration,0);
    sq.ddsTrigDelay = timeAtDrop + TOF - triggerDelay;
    

    sq.anchor(timeAtDrop);
    MakePulseSequence_Rhys(sq.dds,'k',k,'t0',TOF,'T',1e-3,'width',PulseWidth,'dt',dt,...
        'phase',[0,0,0],'chirp',chirp,'delta',delta,...
        'power1',1*[Ch1_Pratio,0,0],'power2',1*[Ch2_Pratio,0,0],'PulseType','Square');
end


%% Blow away |2,1> and |2,2> atoms that are present from mag trap
sq.anchor(timeAtDrop);
% % % Use imaging beam (much weaker than trapping)
% sq.find('Imaging Amp TTL').set(1);
% sq.find('Imaging Freq').set(RunConversions.imaging(0));
% sq.delay(5e-3);
% sq.find('Imaging Amp TTL').set(0);
% sq.anchor(timeAtDrop);

%% Microwave/Raman Stuff
% % % Inputs
MWDelay = 4e-3;
MWDuration = 536e-6;
MagDelay = 0*1e-3;

SGDelay = 17e-3;
SGDuration = 2e-3;

BlowAway1Delay = 1e-3;
BlowAway1Duration = 2e-3;

% % % % % % % % % % % % 
TwoStateImaging = 0;
BiasPrep = 10e-3;
BiasDelay = 2.5e-3;
EWBias = 5;
NSBias = 0;
% EWBias = 0.6;
% NSBias = 7;

RamanTOF = 16.5*1e-3;
RamanPulseWidth = 50*1e-6;
dt = 10e-6;
BeamPower1 = 1;
BeamPower2 = BeamPower1;
% delta = 0.02 + opt.params;
delta = -5 + opt.params;

triggerDelay = 1e-3;
TriggerDuration = 10e-3;

% % % % % % % % % % % % 

BlowAway2Delay = 0.1e-3;
BlowAway2Duration = 1e-3;

% OnOff
BlowAway1 = 1;
Raman = 1;
BlowAway2 = 1;

if RamanAlignment == 1
    BlowAway1 = 0; Raman = 0; BlowAway2 = 0;
    opt.mw.enable(1) = 0;
    opt.mw.enable_sg = 0;
end

% % Timing error checks
if BlowAway1 + Raman + BlowAway2 ~= 0
    if RamanTOF < (MWDelay + MWDuration) + (BlowAway1Delay + BlowAway1Duration)
        error('Raman pulse occurs before MW/blow away1')
    end
    
    if opt.tof < (RamanTOF + RamanPulseWidth) + (BlowAway2Delay + BlowAway2Duration)
        error('Imaging starts before Raman is complete')
    end
end

% % % Sequence 
% MW from |-1,-1> -> |2,0>
if opt.mw.enable(1) == 1
    sq.anchor(timeAtDrop);
    % Set bias
    sq.find('bias e/w').before(MagDelay,5);
    sq.find('bias u/d').before(MagDelay,0); %%% This could be wrong. Previous operation was at zero VOLTS not amps
    sq.find('bias n/s').set(0);

    % Microwave Transfer
    sq.find('R&S list step trig').set(1);
    sq.delay(MWDelay); % delay to prevent state-changing collisions
    sq.find('state prep ttl').set(1);
    sq.delay(MWDuration);
    sq.find('state prep ttl').set(0);

    % return bias to zero
%     sq.find('bias e/w').after(1e-3,0);
    sq.delay(BlowAway1Delay);
end

% Repump light blow away of |-1,-1> atoms
if BlowAway1 == 1
    sq.find('liquid crystal repump').set(7);
    sq.find('Top Repump Shutter').set(0);
    sq.find('Repump Amp TTL').set(1);
    sq.find('Repump Amp').set(10);
    sq.find('Repump Freq').set(RunConversions.repump_freq(0));
    sq.delay(BlowAway1Duration);

    sq.find('Repump Amp TTL').set(0);
    sq.find('Top Repump Shutter').set(1);
    sq.find('liquid crystal repump').set(-2.22);
end

% Stern-Gerlach pulse to test MW transfer
if opt.mw.enable_sg == 1
    %
    % Apply a Stern-Gerlach pulse to separate states based on magnetic
    % moment.  A ramp is used to ensure that the magnetic states
    % adiabatically follow the magnetic field
    %
    sq.anchor(timeAtDrop + SGDelay);
    sq.find('mot coil ttl').set(1);
    t = linspace(0,SGDuration,20);
    sq.find('3d coils').after(t,convert.mot_coil(sq.linramp(t,0,5.5)));
    sq.find('3d coils').after(t,sq.linramp(t,sq.find('3d coils').values(end),convert.mot_coil(0)));
    sq.delay(2*SGDuration);
    sq.find('mot coil ttl').set(0);
    sq.find('3d coils').set(convert.mot_coil(0));
end

% Raman transfer from |2,0> to |1,0> 
if Raman == 1
    %inputs


    sq.anchor(timeAtDrop + RamanTOF);
    sq.find('bias e/w').before(BiasPrep,EWBias);
    sq.find('Bias N/S').before(BiasPrep,NSBias);

    sq.delay(BiasDelay);
    sq.find('bias e/w').set(0);
    sq.find('Bias N/S').set(0);    

    sq.anchor(timeAtDrop + RamanTOF);
    sq.find('Raman DDS Trig').before(TriggerDuration + triggerDelay,1);
    sq.find('Raman DDS Trig').after(TriggerDuration,0);
    sq.ddsTrigDelay = timeAtDrop + RamanTOF - triggerDelay;

    sq.anchor(timeAtDrop);
    chirp = 25.106258428e6;
    k = 22.731334388721734;
    MakePulseSequence_Rhys(sq.dds,'k',k,'t0',RamanTOF,'T',1e-3,'width',RamanPulseWidth,'dt',dt,...
        'phase',[0,0,0],'chirp',chirp,'delta',delta,...
        'power1',[BeamPower1,0,0],'power2',[BeamPower2,0,0],'PulseType','Square');

    sq.anchor(timeAtDrop + RamanTOF + RamanPulseWidth);
%     sq.find('bias e/w').set(0);

end

% Trapping light blow away |2,0> atoms 
if BlowAway2 == 1
    sq.anchor(timeAtDrop + RamanTOF + RamanPulseWidth + BlowAway2Delay);
    sq.find('3D MOT Amp').set(5);
    sq.find('3D MOT Amp TTL').set(1);
    sq.find('3D MOT Freq').set(RunConversions.mot_freq(0));
    sq.delay(BlowAway2Duration);
    sq.find('3D MOT Amp TTL').set(0);
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
    if TwoStateImaging == 1
        makeRamanImagingSequence(sq,'type',Abs_Analysis_parameters.camera,...
            'tof1',opt.tof,'repump delay',0.5e-3,'tof2',opt.tof+opt.misc.tof2,'cycle time',150e-3,...
            'repump Time',200e-6,'pulse Delay',10e-6,'pulse time',10e-6,...
            'imaging freq',imageVoltage,'image freq2',imageVoltage,'repump freq',4.3,'includeDarkImage',true);
    else
        makeImagingSequence(sq,'type',Abs_Analysis_parameters.camera,'tof',opt.tof,...
            'repump Time',100e-6,'pulse Delay',10e-6,'pulse time',[],...
            'imaging freq',imageVoltage,'repump delay',10e-6,'repump freq',4.2,...
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