function varargout = makeSequence_Feedback_Jordan(varargin)   
%%% Note: 
    % This code has been updated to be similar to makeSequenceRyan for NDI and FBC. 

%% Define Jordan's function/notation:
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% %%% %%% Commands:

myOn  = 1;
myOff = 0;
myAlwaysOn = 1;% Used for MOT load and Abs Imag.

%% Parse input arguments
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Always on:

opt = SequenceOptions('load_time',15,'detuning',0,'tof',20e-3,'redpower',2,...
    'keopsys',2);

if nargin == 1
    if ~isa(varargin{1},'SequenceOptions')
        error('If using only one argument it must of type SequenceOptions');
    end
    opt.replace(varargin{1});
elseif mod(nargin,2) == 0
    opt.set(varargin{:});
elseif mod(nargin - 1,2) == 0 && isa(varargin{1},'SequenceOptions')
    opt.replace(varargin{1});
    opt.set(varargin{2:end});
else
    error('Either supply a single SequenceOptions argument, or supply a set of name/value pairs, or supply a SequenceOptions argument followed by name/value pairs');
end



%% Define imaging params:
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Define the params for the incidient field:

% ImageFreq = FtoV('image',opt.detuning + 0*3.1); %at 3 V imaging field, 3.1 MHz detuning, 2.73 MHz at 10 V
ImageFreq = opt.detuning*0.6238/6 + 8.3;
% ImageFreq = opt.params(1);
dipole_field = 5;
% dipole_field = opt.params(1);
imaging_field = 0.5; %%% Smaller as we want NDI
ImageAmp = 7;

%% Initialize sequence
%%% %%% %%% %%% %%%%%% %%% %%% %%% %%%%%% %%% %%% %%% %%%
%%% Define the sequence: 

sq = initSequence;  %load default values (OLD MOT values are default) 
sq.find('87 imag freq').set(ImageFreq);
% sq.find('87 imag freq').set(8.35); % RH updated this. now it does not
% depend on freq? 
sq.find('87 imag amp').set(8);
% sq.find('Variable wave plate').set(-5); %Set for ND imaging

%% MOT loading
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Notes: 
    % * We use CD channel 0b10 = 2 for loading the MOT. 

%%% Code:
if(myAlwaysOn) % myAlwaysOn
    %%% Turn on light for MOT: 
    sq.find('2DMOT').set(1);
    sq.find('3DMOT').set(1);
    sq.find('87 push').set(1);

    %%% 3D MOT beam settings
    sq.find('3DMOT Freq').set(FtoV('trap',24)); % 26 previous
    sq.find('3DMOT amp').set(TrapPtoV('trap',1));
    
    %%% 3D repump beam settings
    sq.find('87 repump').set(1);
    sq.find('Repump Switch').set(0);
    sq.find('87 repump freq').set(FtoV('repump',0));
    sq.find('87 repump amp').set(TrapPtoV('repump',0.93)); % 1 previous
    
    %%% 3D coil settings
    sq.find('H-Bridge Quad').set(1);
    sq.find('CD bit 0').set(0); 
    sq.find('CD bit 1').set(1);
    sq.find('CD2').set(dBtoV('normal',14)); %Coarse control of 3D coils %11 previous
    sq.find('CD Fine/Fast').set(dBtoV('fine',8)); % fine control of 3D coils
    
    %%% Delay for the load time
    sq.delay(opt.load_time);
    
    %%% Turn off the 2D MOT and coils as well as the push beam
    sq.find('2D MOT Coils').before(10e-3,0);
    sq.find('2DMOT').before(10e-3,0);
    sq.find('87 push').before(10e-3,0);
    sq.find('85 push').before(10e-3,0);
end

%% CMOT sequence
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Notes: 
    % * Apply a compressed MOT sequence to temporarily increase the density 
    % by reducing spontaneous emission.  
    % * We switch to CD channel 0b00 = 0 because it is the fast channel.

%%% Code: 
if(myOn) % myOn
    %%% Define time vector for CMOT:
    Tcmot = 15e-3; %20e-3;
    t = 0:1e-3:Tcmot;
    
    %%% Turn on 3D Coils
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(dBtoV('normal',0));
    sq.find('CD Fine/Fast').set(dBtoV('fine',12)); % 9 previous
    
    %%% Turn on trapping light
    sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',46))); % previous 60
    sq.find('3DMOT amp').set(TrapPtoV('trap',1));
    
    %%% Turn on repump:
    sq.find('87 repump freq').set(FtoV('repump',-8.5)); % previous: -7
    sq.find('87 repump amp').set(TrapPtoV('repump',0.05)); % previous: 0.1
    
    %%% Delay for duration of CMOT:
    sq.delay(Tcmot);
end

%% PGC sequence
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Notes: 
    % * Apply polarization gradient cooling to reduce the temperature of 
    % the atoms.  
    % * We use CD channel 0b00 = 0 as it is the fast-switching channel

%%% Code: 
if(myOn) % myOn
    %%% Define time vector for PGC:
    Tpgc = 15e-3; % previous: 25e-3 RH 
    t = 0:0.1e-3:Tpgc;
    
    %%% Change MF and OF params:
    sq.find('CD fine/fast').set(dBtoV('fine',7));
    sq.find('CD0 Fast').set(dBtoV('normal',0));
    sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',72))); % previous: 70
    sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),TrapPtoV('trap',1.3))); % previous: 1.0
    sq.find('87 repump freq').set(FtoV('repump',-4.8)); % previous: -7.25
    sq.find('87 repump amp').set(TrapPtoV('repump',0.004)); % previous: 0.05
    
    %%% Delay for duration of PGC:
    sq.delay(Tpgc);
end

%% Optical pump atoms into the F = 1 manifold
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Notes:
    % * Turn off repump field so that atoms are optically pumped into the 
    % F = 1 % manifold.

%%% Code: 
if(myOn) % myOn
    %%% Define time vector: 
    Tdepump = 1e-3;
    
    %%% Set the depump parmas to turn the light off in the MOT?
    sq.find('repump switch').set(1); %fiber switch off (it's inverted)
    sq.find('87 repump').set(0);
    sq.find('87 repump amp').set(0);
    sq.find('85 repump').set(0);
    % sq.find('85 repump amp').set(0);
    
    %%% Delay for depump process
    sq.delay(Tdepump);
    sq.find('3DMOT').set(0);
end

%% Load into magnetic trap
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Notes: 
    % * Load into the magnetic trap at a high gradient.  We switch quickly 
    % to a low value and then ramp up to the target value

%%% Code: 
if(myOn) % myOn
    %%% Define time vector for loading into MT:
    Tmagload = 150e-3;
    t = 0:5e-3:Tmagload; % turned off. 

    %%% Define magnetic field params: 
    dBLoad = 110;

    %%% Set the magnetic field params:
    sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
    sq.find('CD Fine/Fast').set(dBtoV('fine',0));
    %  sq.delay(Tmagload);
    
    %%% Define time vector for loading into ODT:
    Toptload = 400e-3;
    t = 0:10e-3:Toptload; % turned off.

    %%% Set the ODT params:
    sq.find('Keopsys MO').set(3.9);
    sq.find('Keopsys FA').after(t,sq.minjerk(t,0,DipolePtoV('Keopsys',5))); % previous: 6 RH ; 5 RT
    sq.find('Redpower TTL').set(1);
    sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',15)));
    sq.find('MOT bias').set(1);
    sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));
    
    %%% Delay for duration of load (largest of two processes):
    sq.delay(max(Toptload,Tmagload));
end

%% RF evaporation
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * Remove hot atoms from the sample using RF transitions between the 
    % trapped |F = 1, m_F = -1> state and the untrapped |F = 1, m_F = 0> 
    % state.  
    % * All frequencies are in MHz

if(myOn) % myOn
    %%% Define RF params: 
    rf_start = 20;
    rf_end   = 1;
    rf_rate  = 2.5;  %MHz/s
    rf_ramp_type = 'lin';
    rf_exp_time_constant = 2;
    
    %%% Define RF time vector from params:
    Tevap = (rf_start - rf_end)/rf_rate;
    t = linspace(0,Tevap,50); 

    %%% Set the RF ramping params:
    sq.find('RF atten').set(1);
    if strcmpi(rf_ramp_type,'exp')
        sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
    elseif strcmpi(rf_ramp_type,'lin')
        sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
    end
    
    %%% Delay for duration of evap
    sq.delay(Tevap);
    
    %%% Turn off RF field:
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
end

%% Take dummy images for NDI
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * The cameras first image is dodgy, so we take this now. 
    % * There is also a reference image that is taken to remove the impact
    % of artefacts on the camera. 

%%% Code
if(myOn) % myOn
    if opt.nd.ref_images > 0
        %%% Define IRL time: 
        time_at_evap_end = sq.time;
        sq.anchor(sq.time - 3);
        sq.camDelay = sq.time - 2;

        %%% Set NDI params:
        makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',5e-6,'cycle time',100e-3,...
            'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_amp,'species',85,'num_images',opt.nd.ref_images,...
            'pulse delay',opt.nd.pulse_delay);
        sq.anchor(time_at_evap_end);
    end
end

%% Ramp down trap
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * We need to ramp the coils down instead of turning them off to avoid
    % both eddy currents, and to avoid damaging the coils in the long term

%%% Code
if(myOn) % myOn
    %%% Define time vector for ramping down coils 
    Trampcoils = 0.9; % previous: 0.5
    t = linspace(0,Trampcoils,51);
    
    %%% Define the ramp-down params:
    dB_weak = 0;
    
    %%% Set the ramp-down params: 
    sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));
    
    %%% Delay for duration of ramp-down process:
    sq.delay(Trampcoils);
    
    %%% Old code (indented):
                % sq.find('CD0 Fast').set(0);
                % sq.find('CD Fine/Fast').set(0);
                % % sq.delay(20e-3 - opt.tof);
end

%% Optical evaporation:
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * Main area that I will change. 

%%% Code: 
if(myOn) % myOn
    %%% Define a time vector for optical evap: 
    Tevap = 4; % 3
    t = linspace(0,Tevap,150);

    %%% Define the ODT params from 'opt' 
    final_dipole.RP = opt.redpower;  % defines how low the trap depth is
    final_dipole.FA = opt.keopsys;
     
    %%% Set ODT params:
    sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),0.48)); % previous: 0.40
    sq.find('Keopsys FA').after(t,sq.expramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',final_dipole.FA),0.39)); % previous: 0.4
    
    %%% Delay for duration of optical evap:
    sq.delay(Tevap);
    
    %%% Old code (indented):
                % T = 200e-3;
                % t = linspace(0,T,51);
                % sq.find('Keopsys FA').after(t,sq.linramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',0.8)));
                % sq.find('RedPower CW').after(t,sq.linramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',1.34)));
                % sq.delay(T);
end

%% Adiabatic ramping: 
if(myOff) % myOn %opt.param3
    %%% Define a time vector for optical evap: 
    TAdiabaticRamp = 200e-3; % 200ms
    t = linspace(0,TAdiabaticRamp,150);

    %%% Define the ODT params from 'opt' 
    % final_dipole_ramp.RP = opt.param1;
    % final_dipole_ramp.FA = opt.param2;

    final_dipole_ramp.RP = 1.5;
    final_dipole_ramp.FA = 1.5;
    
    %%% Set ODT params:
    sq.find('RedPower CW').after(t,sq.linramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole_ramp.RP))); 
    sq.find('Keopsys FA').after(t,sq.linramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',final_dipole_ramp.FA))); 
    
        %%% Delay for duration of ramping:
    sq.delay(TAdiabaticRamp);
end

%%% Rethermalisation time was here, moved to below. 

%% Non-destructive imaging
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 

%%% Code: 
if(myOff)
    %%% Set the NDI params from 'opt':
%     makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_delay,'cycle time',opt.nd.cycle_time,...
%         'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_amp,'species',85,'num_images',opt.nd.num_images,...
%         'pulse delay',opt.nd.pulse_delay);
    makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_delay,'cycle time',opt.nd.cycle_time,...
         'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_amp,'species',85,'num_images',opt.nd.num_images,...
         'pulse delay',opt.nd.pulse_delay);
end

% sq.delay(0.1); % hold in the ODT for 100ms. 
% sq.delay(opt.params); % hold in the ODT for variable hold time (USED in callback)

% %% Hold in trap
% if(myOff)
%     % THold = opt.param2;
%     % sq.delay(THold);
% end

%% thermalisation 
if(myOff)
    TThermalisation = 250e-3; % Speilman had 400ms 
    sq.delay(TThermalisation);
end

%% Drop atoms
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * To drop the atoms turn off the ODT's OF.

%%% Code: 
if(myAlwaysOn)
    %%% Set the current time IRL: 
    timeAtDrop = sq.time;
    
    %%%% Turn the fields off. 
    sq.find('2D MOT Coils').set(0);
    sq.find('3DMOT').set(0);
    sq.find('87 repump amp').set(0);
    sq.find('CD0 Fast').set(0);
    sq.find('CD2').set(0);
    sq.find('CD Fine/Fast').set(0);
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('Redpower CW').set(0);
    sq.find('Redpower TTL').after(100e-6,0);
    sq.find('Keopsys FA').set(0);
    sq.find('Keopsys MO').after(100e-6,0);
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
    % sq.find('Variable wave plate').set(-3.7);   %This value switches to absorption imaging
end

%% Take Absorption Image: 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * Ahh

%%% Code: 
if(myAlwaysOn)
    %%% Record time
    sq.anchor(timeAtDrop);

    %%% If there has been no reference images take them and then drop atoms
    if opt.nd.ref_images == 0
        sq.camDelay = timeAtDrop - 2;
    end
    
    %%% Set AI params:
    makeImagingSequence(sq,'tof',opt.tof,'pulse time',40e-6,'repump delay',100e-6,... % previous pulse time: 100e-6 RH ; 40e-6 RT
        'repump time',100e-6,'cam time',5e-6,'cycle time',100e-3,... % previous cycle time: 40e-3?? RH ; 100e-3 RT
        'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
        'fibre switch delay',1e-3,'imaging_field',dipole_field,'image type','horizontal');
    
    %%% Turn of the OF from the ODT: 
    sq.find('Redpower CW').set(0);
    sq.find('Redpower TTL').after(100e-6,0);
    sq.find('Keopsys FA').set(0);
    sq.find('Keopsys MO').after(100e-6,0);
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
    sq.find('3DMOT').set(0);
    sq.find('87 repump amp').set(0);
    
    %%% Old code (indented): 
                % sq.waitFromLatest(60e-3);
                % makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',5e-6,'cycle time',60e-3,...
                %     'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_amp,'species',85,'num_images',1,...
                %     'pulse delay',10e-6);
end

%% Make sure everything is off: 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * Make sure all of the devices are turned off: 

%%% Code: 
if(myAlwaysOn)
    sq.waitFromLatest(0.25);
    setSafeValues(sq);
end 

%% Compile results: 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Notes: 
    % * 

%%% Code: 
if(myAlwaysOn)
    %%% Compile: 
    if nargout == 0
        r = RemoteControl;
        r.upload(sq.compile);
        r.run;
    else
        varargout{1} = sq;
    end
end

%%% End that closes function
end