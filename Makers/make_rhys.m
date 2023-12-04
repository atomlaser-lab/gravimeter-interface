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
% Goal: Collect as many atoms as possible

sq.find('Imaging Freq').set(convert.imaging(opt.detuning));

% % % 2D MOT Light
sq.find('2D MOT Amp TTL').set(1);
sq.find('Push Amp TTL').set(1);

sq.find('2D MOT Freq').set(7.75);
sq.find('Push Freq').set(9.5);

% % % 3D MOT Light
sq.find('3D MOT Amp TTL').set(1);
sq.find('Repump Amp TTL').set(1);
sq.find('MOT Coil TTL').set(1);
% Trapping
sq.find('3D MOT Freq').set(convert.mot_freq(-15));
sq.find('3D MOT Amp').set(5);
sq.find('Liquid crystal Bragg').set(3);
% Repump
sq.find('Repump Freq').set(convert.repump_freq(-1.5));
sq.find('Repump Amp').set(9);
sq.find('Liquid Crystal Repump').set(7);

% % % Mag Field
sq.find('Bias E/W').set(0);
sq.find('Bias N/S').set(0);
sq.find('Bias U/D').set(0.0667);
sq.find('3D Coils').set(convert.mot_coil(1.6));

% Duration
sq.delay(opt.MOT_LoadTime);


%% CMOT
if opt.CMOT_status == 1
    % Goal: Increase confinement with minimal loss
    % % % Stop blowing hot atoms in
    sq.find('2D MOT Amp TTL').before(10e-3,0);
    sq.find('push amp ttl').before(10e-3,0);

    % % % Reduce Re-radiation Pressure
    sq.find('3D MOT freq').set(convert.mot_freq(-25)); % 26
    sq.find('3D MOT Amp').set(4.5);

    sq.find('repump freq').set(convert.repump_freq(-5)); %-8
    sq.find('Repump Amp').set(9); %9

    % % % Increase Mag Field
    sq.find('3D Coils').set(convert.mot_coil(convert.mot_coil_reverse(sq.find('3D Coils').values(end)) - 0.3));
    sq.find('Bias E/W').set(0);
    sq.find('Bias N/S').set(0);
    sq.find('Bias U/D').set(0.0667);

    % % %Duration
    sq.delay(10e-3); 
end


%% PGC
if opt.PGC_status == 1
    % Goal: Produce cloud cold (with minimal atom loss) enough to be loaded into mag trap
    PGC_Time = sq.time;
    % Durations
    tStep = 0.5e-3;
    t_TrappingRamp = 40e-3;
    t_RepumpRamp = 5e-3;
    t_MagRamp = opt.params*1e-3;
    t_MagDelay = 0*1e-3;

    % % % Set magnetic field to zero
    sq.anchor(PGC_Time - t_MagDelay);
    t = linspace(0,t_MagRamp,100);
    sq.find('Bias E/W').after(t,sq.linramp(t,sq.find('Bias E/W').values(end),0));
    sq.find('Bias N/S').after(t,sq.linramp(t,sq.find('Bias N/S').values(end),0));
    sq.find('Bias U/D').after(t,sq.linramp(t,sq.find('Bias U/D').values(end),2.707));
    sq.find('3D Coils').after(t,sq.linramp(t,sq.find('3D Coils').values(end),convert.mot_coil(0)));
    sq.delay(t_MagRamp);

    % % % Reduce heating 
    sq.anchor(PGC_Time);
    t = linspace(0,t_RepumpRamp,round(t_RepumpRamp/tStep,0));
    sq.find('Repump Amp').after(t,sq.linramp(t,sq.find('Repump Amp').values(end),7));
    sq.find('repump freq').after(t,sq.linramp(t,sq.find('repump freq').values(end),convert.repump_freq(-8.75))); 
    sq.delay(t_RepumpRamp);
 
    sq.anchor(PGC_Time);
    t = linspace(0,t_TrappingRamp,round(t_TrappingRamp/tStep,0));
    sq.find('3D MOT Amp').after(t,sq.linramp(t,sq.find('3D MOT Amp').values(end),3.9));
    sq.find('3D MOT Freq').after(t,sq.linramp(t,sq.find('3D MOT Freq').values(end),convert.mot_freq(-70)));
    sq.delay(t_TrappingRamp);
end

%% Mag Load
if opt.LoadMagTrap_status == 1
    % % % Pump into |F=1, mf = -1>
    % Repump off
    sq.find('Repump Amp TTL').set(0);
    sq.find('Top repump shutter').set(1);
    sq.find('liquid crystal repump').set(-2.22);
    % Trapping amplitude
    sq.find('3D MOT Freq').set(convert.mot_freq(-70)); %-71
    sq.find('3D MOT Amp').set(4.5); %3.9

    sq.delay(2e-3);

    % % % Ramp on Mag Field
    

    % % % Turn Light Off


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