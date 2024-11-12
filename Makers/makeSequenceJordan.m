function varargout = makeSequenceJordan(varargin)   

%% Prelude: (No code):
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Prelude:
    % * Welcome to my sequence. This code defines the sequence to make a BEC, and provides comments on each stage, 
    % and element. The following notations are used. 
    % * Notation: 
        % '*'  -- A dot point is used to signify a single bit of information, which may carry over several lines. 
        % '->' -- Is a particular element, and its purpose
        % '%%' -- Indicates a new section, and forces MATLAB to place a vertical line in the code. 
        % 'Comment' -- A comment on previously used values, or nominal values of a paramter. 
        % 'xx' -- Placeholder name
    % * All elements are triggered by a voltage, however we want to know the mapping to a physical parameter, so we 
    % define function: 'XtoV('__Name_of_element__', __Value_of_element__), which maps the parameter symbolised by 
    % X to a voltage control. 

%% Parse input arguments (ALWAYS ON):
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Parse input arguments:
    % This will define the opt parameters -- when the code is opened these are the default values that are used.
opt = SequenceOptions('load_time', 5, 'detuning', 0, 'tof', 20e-3, 'redpower', 2, 'keopsys',2);
    % Previously known values of the laser powers to make a BEC:
        % keyopsis bec 0.8 0.75
        % redpower bec 1.34

    % %%% %%% I AM NOT SURE WHAT THIS CODE DOES:
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

    % Define the imaging parameters:  
ImageFreq = opt.detuning*0.6238/6 + 8.3; % ??
% ImageFreq = opt.params*0.6238/6 + 8.3; % Comment: Scan detuning over range set in params
% ImageFreq = 0;
dipole_field = 5; % ??
% imaging_field = 0.5;
imaging_field = 9.25; % Previous value: 5
ImageAmp = 5; % Previous value: 5 -- until and including mag load

%% 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Initialize sequence:
    % * All elements are controlled with voltages that have names defined in sq, we use sq.find('__Name_of_element__') 
    % to find the element, and then apply .set(__Value_of_element__) to specify a value. The full command has the 
    % form: sq.find('__Name__').set(__Value__). Commands seen below:

sq = initSequence; % Command loads the default values for the OLD MOT values. 
sq.find('87 imag freq').set(8.35); % Defines the imaging lights frequency.
sq.find('87 imag amp').set(8.00);  % Defines the imaging lights amplitude.

%% My Commands (No code)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% My commands:
    % Define commands for turning sections on and off, as I have a small brain. This avoids mistakes. 
myOn = 1;
myOff = 0;
    % Certain parts of the sequence should always be on, so we define a parameter to be always on. 
myOnAlways = 1;

%% MOT Loading Sequence (ALWAYS ON):
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% MOT loading sequence:
    % Comment: We use CD channel 0b10 = 2 for loading the MOT.
    % * When making a MOT we need the following things: mag field, and light.
    % * As a flow chart this is how our system should work with parameters listed next to each stage:
        % -> Load time: The load time is specified in the 'opt' command,  and is changed in the command 
        % line via, opt.load_time = __. XXXX__ This will later be changed to time_load, to abide by the latter
        % used notation. 
        % -> Time of MOT stage: time_MOT -- how long should the MOT stage be applied for?
        % -> Dispenser on (always) -- no need to worry about this. 
        % -> 2DMOT_light: nu_2D, A_2D -- now we need atoms to go to the 3DMOT, this is down with push beam.
        % -> 2DMOT_mag: Current driver (CD) and H-bridge need to be triggered.
        % -> Push beam: nu_push, A_push.
        % -> 3DMOT_light: nu_3D, A_3D.
        % -> 3DMOT_mag: Current driver (CD) and H-bridge need to be triggered.
        % -> Repump: nu_repump, A_repump. 
        % -> Once the atoms have loaded into the 3DMOT from the 2DMOT I need to wait a finite amount of time, 
        % 'wait_time_2D' before turning the 2DMOT, Push beam and other things off to make sure optimal transfer
        % between 2DMOT and 3DMOT.
    % * After this stage assuming a 15 s load time, we should see __ atoms with a temperature of __ uK after 
    % release from the MOT sequence.  
    % COMMENTS:
        % 2DMOT Power, frequencey, and magnetic field are constant 
    % * Now we define the parameters for the MOT loading sequence:

%%% Sequence:
if(myOnAlways) % myOn myOff
    %%% Turn on the 2D, 3D, and push beam
    sq.find('2DMOT').set(1);
    sq.find('3DMOT').set(1);
    sq.find('87 push').set(1);

    %%% Trapping light:
    sq.find('3DMOT Freq').set(FtoV('trap',23)); %23 %15.7412
    sq.find('3DMOT amp').set(TrapPtoV('trap',1));

    %%% Repump beam: 
    sq.find('87 repump').set(1);
    sq.find('Repump Switch').set(0);
    sq.find('87 repump freq').set(FtoV('repump',0)); %0
    sq.find('87 repump amp').set(TrapPtoV('repump',1)); %1

    %%% 3D Magnetic field coils: 
    sq.find('H-Bridge Quad').set(1);
    sq.find('CD bit 0').set(0); 
    sq.find('CD bit 1').set(1);
    sq.find('CD2').set(dBtoV('normal',11)); %Coarse control of 3D coils %11
    sq.find('CD Fine/Fast').set(dBtoV('fine',8)); % fine control of 3D coils

    %%% Delay for the load time
    sq.delay(opt.load_time);

    %%% Turn off the 2D MOT and coils as well as the push beam
    sq.find('2D MOT Coils').before(10e-3,0);
    sq.find('2DMOT').before(10e-3,0);
    sq.find('87 push').before(10e-3,0);
    sq.find('85 push').before(10e-3,0);

    %%% Additional command -- probably delete this at some point. 
    sq.find('C6 - N/C').set(1).after(100e-6,0); %setup SLM -- PROBS DONT NEED THIS
end

%% CMOT Sequence (ON)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% CMOT sequence:
    % Comment: Use CD channel 0b00 = 0, as this is the fast control.
    % * To increase the PSD, we need to increase the density via: \rho_\phi = n \lambda_{dB}^3 which implies 
    % that: V \propto T^{3/2} \propto 1/n.
    % Additionally, the increased density will decrease the spontaneous emission rate. 
    % * Next, when we refer to MOT this is the 3D as the 2D MOT is only used above to load into the 3D MOT.  
    % * Again, we want to make a flow chart: 
        % -> Time of CMOT stage: time_CMOT -- how long should the CMOT stage be applied for?
        % -> 3DMOT_light: {nu_3D, A_3D} --> {nu_3D, A_3D} -- we want to change the lights to change the frequency,
        % and amplitude. 
        % -> 3DMOT_mag: The CD and H-bridge needs to be triggered 
        % -> Repump_light: To make sure our atoms stay in the correct |F, m_F> state we need to use a repump light. Which 
        % needs to have: nu_repump, A_repump. 
        % -> 
    % * After the CMOT stage and assuming we started with __ atoms after a 15 s load time, we should see __ atoms with 
    % a temperature of __ uK after CMOT. 
    % * Now we define the parameters for the CMOT loading sequence:

%%% Sequence:
if(myOn) % myOn myOff
    %%% Define the time vector:
    Tcmot = 15e-3; %15e-3
    t = 0:1e-3:Tcmot;

    %%% 3D Magnetic field coils:
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(dBtoV('normal',0));
    sq.find('CD Fine/Fast').set(dBtoV('fine',9)); 

    %%% Trapping light
    sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',45))); %45
    sq.find('3DMOT amp').set(TrapPtoV('trap',1));

    %%% Repump beam:
    sq.find('87 repump freq').set(FtoV('repump',8)); %-4
    sq.find('87 repump amp').set(TrapPtoV('repump',0.2)); %0.8

    %%% Delay steps for duraction of T_CMOT.
    sq.delay(Tcmot);
end

%% PGC Sequence (ON)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% PGC sequence:
    % Comment: Use CD channel 0b00 = 0, as this is the fast control.
    % * Using PGC decreases the temperature of the atoms.
    % * We interrogate the atoms using the repump light in addition to the
    % light being held in the 3DMOT. 

    % * Now we define a flow chart of the steps involved in PGC:
        % -> Time of PGC stage: time_PGC -- how long should PGC be applied for?  
        % -> 3DMOT_light (for trapping): nu_3D, A_3D
        % -> 3DMOT_mag: The CD and H-bridge needs to be triggered 
        % -> Repump_light (for cooling (PGC)): nu_repump, A_repump
    % * After this stage and assuming we started with __ atoms after a 15 s load time, we should see __ atoms with a 
    % temperature of __ uK after PGC. 
    % * Now we define the parameters for the PGC loading sequence:

%%% Sequence:
if(myOn) % myOff
    %%% Define the time vector:
    Tpgc = 25e-3; %22e-3
    t = 0:1e-3:Tpgc;

    %%% 3D Magnetic field coils:
    sq.find('CD fine/fast').set(dBtoV('fine',7));
    sq.find('CD0 Fast').set(dBtoV('normal',0));

    %%% Trapping light:
    sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',70))); %56 %70
    sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),TrapPtoV('trap',1)));

    %%% Repump beam:
    sq.find('87 repump freq').set(FtoV('repump',-9)); %-2.5 
    sq.find('87 repump amp').set(TrapPtoV('repump',0.01)); %0.05
    
    %%% Delay for the duraction of PGC:
    sq.delay(Tpgc);
end

%% Optical pump into |F = 1> (ON):
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Optical pump atoms into the F = 1 manifold:
    % Comment: 
    % * After each stage the atoms rethermalise and redistribute themselves. This means atoms are in all the manifolds. 
    % We need to return atoms the the |F = 1> manifold.  Turn off repump field so that atoms are optically pumped into 
    % the F = 1 manifold.
    % * For magnetic trapping (next stage) we need atoms in the |1, -1> state, and hence we want as many atoms in the 
    % |F = 1> manifold. 

%%% Sequence:
if(myOn) % myOff
    %%% Define time vector:
    Tdepump = 1e-3;

    %%% Fibre switch??
    sq.find('repump switch').set(1); %fiber switch off (it's inverted)

    %%% Trapping light:
    sq.find('87 repump').set(0);
    sq.find('87 repump amp').set(0);

    %%% HMM WHY ARE WE USING THE 85 LASER?
    sq.find('85 repump').set(0);
    sq.find('85 repump amp').set(0);
    
    %%% Delay for duration of optical pumping:
    sq.delay(Tdepump);

    %%% Turn off the 3D MOT:
    sq.find('3DMOT').set(0);
end

%% Load into Magnetic trap (ON)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Load into magnetic trap:
    % Comment: 
    % Load into the magnetic trap at a high gradient. We switch quickly to
    % a low value and then ramp up to the target value.
    % 

%%% Sequence:
if(myOn) % myOn myOff
    %%% Define a time vector:
    Tmagload = 150e-3; %150e-3
    t = linspace(0,Tmagload,50);

    %%% Define magnetic field gradient (dBLoad)
    dBLoad = 110; %110
    sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
    sq.find('CD Fine/Fast').set(dBtoV('fine',0));
    sq.delay(Tmagload); %Tmagload

    %%% %%% %%% Dipole trap:

    %%% Define a time vector:
    Toptload = 400e-3; %400e-3
    t = linspace(0,Toptload,50);

    %%% Turn on the dipole beams:
    sq.find('Keopsys MO').set(3.9);
    sq.find('Keopsys FA').after(t,sq.minjerk(t,0,DipolePtoV('Keopsys',5))); %5
    sq.find('Redpower TTL').set(1);
    sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',15))); %15
    
    %%% Turn on the magnetic field biases:
    sq.find('MOT bias').set(1);
    sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));
    
    %%% Delay for duration of sequence. 
    sq.delay(max(Toptload,Tmagload));
end

%% RF Evaporation (ON)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% RF evaporation:
    % Comment: 
    % Remove hot atoms from the sample using RF transitions between the trapped |F = 1, m_F = -1> 
    % state and the untrapped |F = 1, m_F = 0> state.  
    % All frequencies are in MHz 
    % Start with a a low value and then ramp up to the target value.

%%% Sequence:
if (myOn) % myOn myOff
    %%% Define the start and end points for the RF frequencies. 
    rf_start = 20; % MHz
    rf_end = 1;    % MHz Note that rf_end == 1 if dipoles are correctly positioned
        
    %%% Define the ramp rate and the type of ramp. 
    rf_rate = 2.2; %3 %MHz/s
    rf_ramp_type = 'lin';
    rf_exp_time_constant = 2;
    
    %%% Define the time vector:
    Tevap = (rf_start - rf_end)/rf_rate;
    t = linspace(0,Tevap,100);

    %%% Compare the rf_ramp_type and do different actions for linear and exponential ramping. 
    sq.find('RF atten').set(1);

        %%% if exponential ramping:
    if strcmpi(rf_ramp_type,'exp') 
        %%% Ramp function: freq_ramp_t = freq_ramp_init * exp(t / tau) --> if t = end(t), then freq_ramp_t = freq_ramp_fin. 
        sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
        %%% if linear ramping:
    elseif strcmpi(rf_ramp_type,'lin')
        %%% Ramp function: freq_ramp_t = freq_ramp_int - freq_ramp_ramp * t --> if t = end(t), then freq_ramp_t = freq_ramp_fin. 
        sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
    end
    
    %%% Delay for duration of sequence:
    sq.delay(Tevap);
    
    %%% Turn off the signal:
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));

    %%% Old comments: Turn off magnetic trap (use for dipole alignment)
    % sq.find('CD0 Fast').set(0);
    % sq.find('CD Fine/Fast').set(0);
    % sq.delay(20e-3 - opt.tof);
end



%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%% Take dummy images | Section may not be needed at the start!!
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Take dummy images:
    % Comment: 
% if opt.nd.ref_images > 0 % Original line used by RT. 
if(myOff) % myOn myOff
    %%% Define the time the images are taken?:
    time_at_evap_end = sq.time;
    sq.anchor(sq.time - 3);
    sq.camDelay = sq.time - 2;

    %%% Make the NDI imaging sequence:
    makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',5e-6,'cycle time',100e-3,...
        'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_amp,'species',85,'num_images',opt.nd.ref_images,...
        'pulse delay',opt.nd.pulse_delay);
    sq.anchor(time_at_evap_end);
end
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%



%% Turn off coils (ON):
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Ramp down coils:
    % Comment: 
    % * This section is to slowly turn off the coils. 
    % * Most importantly, we want as many atoms transferred from the  magnetic trap to the dipole trap. By slowly 
    % decreasing the size of  the magnetic trap we increase the chances  that atoms fall into the dipole trap. 
    % * Doing this slowly will result in less eddy currents forming. Due to d_t B we have a spatially varying 
    % electric field (-curl(E)). This can temporally vary, and induce a spatially changing magnetic field via 
    % d_t E = c^2(curl(B) - mu_0 j). Hence causing changes in magnetic field. Slowly ramping down the coils 
    % increases dt so d_t B ~ 0, and same for d_t E ~ 0. 

%%% Sequence:
if(myOn) % myOn myOff
    %%% Define time vector:
    Trampcoils = 0.5;
    t = linspace(0,Trampcoils,51);
    
    %%% We ramp down from dB_current to dB_weak == 0:
    dB_weak = 0;
    sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));

    %%% Delay for duration of sequence:
    sq.delay(Trampcoils);
    
    %%% Old commands: 
    % rf_final = 1;
    % sq.find('RF frequency').after(t,sq.linramp(t,sq.find('Rf frequency').values(end),RFtoV(rf_final)));
    % sq.delay(1); %to hold in dipole trap
end



%%% (UNSURE) 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%% Blow away |F = 2> | Section may not be needed at the start!!
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Blow away |F = 2> manifold:
    % Comment: 
    % * Atoms will again be an a distribution of states -- we will remove
    % atoms in the |F = 2> manifold. 
    % * WHYYYY do we use the 87 imaging light?? 

%%% Sequence:
if(myOff) % myOn myOff
    %%% Turn on imaging light to push away |F = 2> manifold, flash, and then turn off. 
    sq.find('87 imag').set(1);
    sq.delay(1e-3);
    sq.find('87 imag').set(0);
end
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% (UNSURE) 



%%% (UNSURE) 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%% HH configuration | Section may not be needed at the start!!
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Switch to Helmholtz configuration for state preparation:
    % Comment:
    % * Figure out why we need a Helmholtz field configuration for
    % preparing atoms? 

%%% Sequence:
if(myOff) % myOn myOff
    %%% Turn off magnetic trap entirely 
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);

    %%% Ensure that the quadrupole (anit-Helmholtz) configuration is off.
    sq.find('H-Bridge Quad').set(0);

    %%% Wait time to ensure that both configurations are not simultaneously triggered. 
    sq.delay(50e-6);

    %%% Turn on Helmholtz configuration. 
    sq.find('H-Bridge Helm').set(1);

    %%% Delay for duration of sequence:
    sq.delay(100e-6);

    %%% Define a time vector for ramp
    Tramp = 1;
    t = linspace(0,Tramp,51);

    %%% Turn on a small magnetic field gradient to prepare state?? 
    % WHYYYYYYYYYYYY?
    sq.find('CD0 Fast').after(t,sq.minjerk(t,0,dBtoV('normal',20)));
    %sq.find('MOT bias coil').after(t,sq.minjerk(t,sq.find('MOT bias coil').values(end),0));
    %sq.delay(Tramp);
end
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% (UNSURE) 



%% Optical evaporation (ON): 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Optical evaporation:
    % Comment:
    % * in RT's code he uses Tevap = 2; however I am evaporating for longer. I am not sure why?? 
    % * RT uses 0.25 instead if 0.42 ... what does this control? 

if(myOn) % myOn myOff
    %%% Define a time vector:
    Tevap = 4; % 3
    t = linspace(0,Tevap,150);

    %%% Define the dipole powers in terms of the opt params:
    final_dipole.RP = opt.redpower;
    final_dipole.FA = opt.keopsys; 

    %%% Exponentially ramp up the dipole beam powers to opt params values:
    sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),0.42)); %0.42
    sq.find('Keopsys FA').after(t,sq.expramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',final_dipole.FA),0.42)); %0.42

    %%% ??? Why do we need a VWP? 
    sq.find('variable wave plate').set(-2.42);

    %%% Delay for duration of evaporation:
    sq.delay(Tevap);
    
    %%% ??? Why do we delay all by 1?? 
    sq.delay(1);
end

%% NDI (OFF)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Non-destructive imaging (NDI):
    % Comment:
    % * TBH I am not sure what is really going on here, so I will need to figure this section out soon. 

%%% Sequence:
if(myOff) % myOn myOff
    makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_delay,'cycle time',opt.nd.cycle_time,...
        'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_amp,'species',85,'num_images',opt.nd.num_images,...
        'pulse delay',opt.nd.pulse_delay);
end

%% Drop atoms (ALWAYS ON):
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Drop atoms 
%%% (ALWAYS ON) 
    % Comments:
    % * This stage we want to know the real time (since the start), and then we want to turn all 
    % of the electronics off (excluding imagining). 

% sq.delay(500e-3);

%%% Sequence:
if(myOnAlways)
    %%% Determine the drop time:
    timeAtDrop = sq.time;

    %%% Turn of EVERYTHING (besides imagining):

        %%% Magnetic fields:
    sq.find('2D MOT Coils').set(0);
    sq.find('3DMOT').set(0);
    sq.find('CD0 Fast').set(0);
    sq.find('CD2').set(0);
    sq.find('CD Fine/Fast').set(0);
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    
        %%% Optical fields:
    sq.find('87 repump amp').set(0);
    sq.find('Redpower CW').set(0);
    sq.find('Redpower TTL').after(100e-6,0);
    sq.find('Keopsys FA').set(0);
    sq.find('Keopsys MO').after(100e-6,0);

        %%% RF fields:
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
% sq.find('Variable wave plate').set(-3.7);   %This value switches to absorption imaging
end 

%% AbsImag (ALWAYS ON)
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Take Absorption Image
%%% (ALWAYS ON) 
    % Comments:
    % * This should always be on

%%% Sequence:
if(myOnAlways)
    %%% Link up time from when atoms are dropped
    sq.anchor(timeAtDrop);
    if opt.nd.ref_images == 0
        sq.camDelay = timeAtDrop - 2;
    end

    makeImagingSequence(sq,'tof',opt.tof,'pulse time',100e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',50e-6,'cycle time',100e-3,...
    'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
    'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal');

    %%% Turn off the dipoles
    sq.find('Redpower CW').set(0);
    sq.find('Redpower TTL').after(100e-6,0);
    sq.find('Keopsys FA').set(0);
    sq.find('Keopsys MO').after(100e-6,0);
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
    sq.find('3DMOT').set(0);
    sq.find('87 repump amp').set(0);

    % sq.waitFromLatest(60e-3);
    % makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',5e-6,'cycle time',60e-3,...
    %     'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_amp,'species',85,'num_images',1,...
    %     'pulse delay',10e-6);
end

%% 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Wait sequence to return electronics to safe values:
sq.waitFromLatest(0.25);
setSafeValues(sq);

%% Ending sequence (ALWAYS ON):
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% End sequence
    % Comments:
    % * This just compiles, uploads and completes the code, before the function is ended from line 1. 

%%% Sequence:
if(myOnAlways)
    if nargout == 0
        r = RemoteControl;
        r.upload(sq.compile);
        r.run;
    else
        varargout{1} = sq;
    end
end

%%% End the function (from line 1):
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
    %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
        %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
            %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
                %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
                    %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
                        %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
                    %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
                %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
            %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
        %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
    %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 