function varargout = makeSequence_Feedback_Jordan_NEW2(varargin)

%%% --- --- --- %%% Preamble to set-up the code:
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

%%
%%% --- --- --- %%% Set-up stages:
    %%%  Camera Type
Camera = "high res";

%%% --- --- --- %%% Cooling sequence
    %%% Not if stage _j_ is off, all stages after _j_ are off. 
CMOT = 1  ;
PGC = 1  ;
Load_Mag_Evap = 1  ;
RF_Evaporation = 1  ;
Load_Opt_Evap = 1  ;
Optical_Evaporation = 1  ;

%%% --- --- --- %%% Feedback related stages
    %%% These can be toggled on/off independantly. 
Hold_In_ODT = 0  ; %%% only for determing trap life-time.
Ramp_Up_Powers = 1  ;
Shadowgraph_Imaging = 1  ;
Thermalisation_Time = 0  ;

    %%% Display Information Status
Display_info = 1 ;

%%% --- --- --- %%%  Initialisation of useful conversion for dipoles and imaging
    %%% Conversion factor from detuning to frequency of the imaging beam
ImageFreq = (opt.detuning+43.5)^4*6e-8-(opt.detuning+43.5)^3*8e-7-(opt.detuning+43.5)^2*0.0004+(opt.detuning+43.5)*0.0796+5.4703; %assuming imaging field is set to 5V
    
    %%% Field strngth:
dipole_field = 5; %
imaging_field = dipole_field; %5
ImageAmp = 9; %5 until and including mag load, 7 for high OD samples on horizontal imaging

%%
%%% --- --- --- %%% Initialisation of the MOT sequence
    %%% Define the time of the MOT:
Tmot = opt.load_time;

    %%% Initialise the sequence:
sq = initSequence; %load default values

    %%% Turn on the trapping light
sq.find('87 imag freq').set(8.3); % 8.35
sq.find('87 imag amp').set(8); % 8

    %%% Turn on the MOT light and push beam
sq.find('2DMOT').set(1);
sq.find('3DMOT').set(1);
sq.find('87 push').set(1);

    %%% 3D MOT beam settings
sq.find('3DMOT Freq').set(6.6); % 6.7 %*was 6.6*
sq.find('3DMOT amp').set(4.4); %8 %*was 4.4*

    %%% 3D repump beam settings
sq.find('Repump Switch').set(0);
sq.find('87 repump').set(1);
sq.find('87 repump freq').set(4.4); %4.68 %*was 4.4*
sq.find('87 repump amp').set(7); %*was 7*
    
    %%% 3D coil settings
sq.find('H-Bridge Quad').set(1);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(1);
sq.find('CD2').set(1.80); %1.8 %*was 1.8*
sq.find('CD Fine/Fast').set(0); %8 fine control of 3D coils

    %%% Earth canceling biases
sq.find('Earth Bias 1').set(3); %6
sq.find('Earth Bias 2').set(5); %4
sq.find('Earth Bias 3').set(0); %2

    %%% Delay for sequence duration
sq.delay(Tmot);

%%
%%% --- --- --- %%%  Compression MOT sequence
if CMOT == 1
    %%% Turn off the coils, and the push beam. 
sq.find('2D MOT Coils').before(10e-3,0);
sq.find('2DMOT').before(10e-3,0);
sq.find('87 push').before(10e-3,0);
sq.find('85 push').before(10e-3,0);

    %%% Time for CMOT to be applied
Tcmot = 30e-3; %10e-3
t = 0 : 2e-3 : Tcmot;

    %%% 3D Coils 
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('CD0 Fast').set(dBtoV('normal',0));
sq.find('CD Fine/Fast').set(7.5); %2

    %%% Trapping light
% sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',35))); %45
sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),4.6)); %4.5
sq.find('3DMOT amp').set(6); %4.5

    %%% Repump light 
% sq.find('87 repump freq').set(FtoV('repump',6)); %-7.5 %-4
sq.find('87 repump freq').set(6.25); %6.2
sq.find('87 repump amp').set(6); %5

    %%% Delay for sequence duration
sq.delay(Tcmot);

%%
%%% --- --- --- %%% PGC
if PGC == 1
    %%% Time for PGC to be applied:    
Tpgc = 6e-3; %5e-3
t = 0 : 1e-3 : Tpgc;

    %%% Turn on the earth cancelling biases
if(1)
    sq.find('Earth Bias 1').set(1.2);  % Probably N/S
    sq.find('Earth Bias 2').set(6.5);  % Probably U/D
    sq.find('Earth Bias 3').set(0.05); % Probably E/W
end 

    %%% Turn off the other magnetic fields:
sq.find('CD fine/fast').set(0); %0.2
sq.find('CD0 Fast').set(0); %1
    
    %%% Reduce the optical field to decrease radiation pressure from MOT
sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),2));
sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),4.1)); %4

    %%% Reduce the optical field to decrease radiation pressure from 87 trapping
sq.find('87 repump freq').set(2.3); %2.25
sq.find('87 repump amp').set(3); %3

    %%% Delay for sequence duration
sq.delay(Tpgc);
 
%%
%%% --- --- --- %%%  Optical pumping into |F  = 1> manifold
Tdepump = 1e-3;
sq.find('repump switch').set(1); %fiber switch off (it's inverted)
sq.find('87 repump').set(0);
sq.find('87 repump amp').set(0);
% sq.find('85 repump').set(0);
sq.delay(Tdepump);
sq.find('3DMOT').set(0);

%%
%%% --- --- --- %%%  Load into the Magnetic Trap
if Load_Mag_Evap == 1
    %%% Time to load into the magnetic trap:    
Tmagload = 150e-3; % prev 200e-3

    %%% Turn off earth canceling biases:
sq.find('Earth Bias 1').set(0);
sq.find('Earth Bias 2').set(0);
sq.find('Earth Bias 3').set(0);

    %%% Turn on magnetic field
t =  0:5e-3:Tmagload;
dBLoad = 110; % opt 110 
sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',55),dBtoV('normal',dBLoad)));
sq.find('CD Fine/Fast').set(dBtoV('fine',4));
% sq.delay(Tmagload); %Tmagload
    
    %%% Turn on dipoles to load into them
Toptload = 400e-3; % prev 100e-3
t = linspace(0,Toptload,40);
sq.find('Keopsys FA').after(t,sq.minjerk(t,0,DipolePtoV('keopsys',3))); %8  %2.6 (when no PtoV used) %%%GOES TO 3.05 AT MOST % 7.5 W. 
sq.find('Redpower TTL').set(1);
sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',11))); %14 % 15
sq.find('MOT bias').set(1);
sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));

    %%% Delay for duration of sequence
sq.delay(max(Toptload,Tmagload));
 
 
    %%% *Delay in the magtrap to see the life time (for testing purposes only!)*
        % * Using sq.delay(1) holds atoms in the mag trap for 1s
        % * Using sq.delay(opt.params) holds atoms in the mag trap for a variable time. 
        % This means going into the callback "Callback_MeasureTemperature_fancy.mlx" and 
        % changing the variable to opt.params

% sq.delay(0); % delay for n second(s).
% sq.delay(opt.params); %0 to as long as it can holds

%%
%%% --- --- --- %%% RF Evaporation
if RF_Evaporation == 1
    %%% RF evapiration start point, end point, rate, and time
rf_start = 16; % opt 16
rf_end = 1.5; % opt 1.5
rf_rate = 4.5; %MHz/s % opt 4
Tevap = (rf_start - rf_end)/rf_rate;
rf_ramp_type = 'lin';
rf_exp_time_constant = 2;
t = linspace(0,Tevap,25);
    
    %%% Ramp the RF freq down:
sq.find('RF atten').set(1);
if strcmpi(rf_ramp_type,'exp')
    sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
elseif strcmpi(rf_ramp_type,'lin')
    sq.find('RF frequency').set(RFtoV(rf_start)); % *NEW IN RT's SCRIPT *
    sq.delay(0.25);                               % *NEW IN RT's SCRIPT *
    sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
end

    %%% Delay for duration of sequence
sq.delay(Tevap);

    %%% Turn RF off
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));
 
    %%% TESTING: 
        % ** Test for Magnetic trap and optical dipole alignment (to be turned on only 
        % in testing procedures- check within few ms tof).*
% sq.find('CD0 Fast').set(0);
% sq.find('CD Fine/Fast').set(0);
% sq.delay(20e-3 - opt.tof);

        % ** Turn off dipoles before ramping down coils (to be uncommented only for 
        % testing procedure of the atoms in the magtrap)*
        % sq.find('Redpower CW').set(0);
        % sq.find('Redpower TTL').after(100e-6,0);
% sq.find('Keopsys FA').set(0);
% sq.find('Keopsys MO').after(100e-6,0);

%%
%%% --- --- --- %%% Take dummy images:
if Shadowgraph_Imaging == 1
    if opt.nd.ref_images > 0
        %%% Define IRL time: 
        time_at_evap_end = sq.time;
        sq.anchor(sq.time - 3);
        sq.camDelay = sq.time - 2;
        
        %%% Set NDI params:
        makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',5e-6,'cycle time',100e-3,...
            'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_power,'species',85,'num_images',opt.nd.ref_images,...
            'pulse delay',opt.nd.pulse_delay);
        
        %%% Delay for duration of sequence
        sq.anchor(time_at_evap_end);
    end
end 

%% 
%%% --- --- --- %%% Load into Optical Dipoles - Ramp down coils
    % Note: 
        % * For slow ramp down use 0.9 (s).
        % * For fast ramp down use 0.3 (s).

if Load_Opt_Evap == 1
    %%% Ramp the coils down
Trampcoils = 0.8; %0.4 opt
t = linspace(0,Trampcoils,25);
dB_weak = 0;
sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));
% sq.find('CD Fine/Fast').after(t,sq.linramp(t,sq.find('CD Fine/Fast').values(end),0));
    
    %%% Delay for duration of sequence
sq.delay(Trampcoils);

%% 
%%% --- --- --- %%% Optical Evaporation
if Optical_Evaporation == 1
    %%% Time for optical evap    
Tevap = 4; % 2.8 opt

    %%% Set the evaporation end points for the dipole beams:
final_dipole.RP = opt.redpower; % 1.00 opt
final_dipole.FA = opt.keopsys;  % 1.95 opt

    %%% Begin the ramp-down sequence
t = linspace(0,Tevap,100);
sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP), 0.4)); % opt 0.55
sq.find('Keopsys FA' ).after(t,sq.expramp(t,sq.find('Keopsys FA' ).values(end),DipolePtoV('keopsys' ,final_dipole.FA), 0.4)); % opt 0.7   
% sq.find('variable wave plate').set(-2.42);
    
    %%% Delay for duration of sequence
sq.delay(Tevap);

 
%%% --- --- --- %%%  End the cooling sequence:
end %%% Optical evaporation
end %%% Load in optical dipoles
end %%% RF evaporation
end %%% Load into the Mag trap
end %%% PGC
end %%% CMOT

%%% --- --- --- %%% %%% --- --- --- %%% %%% --- --- --- %%% %%% --- --- --- %%% %%% --- --- --- %%% 
%%% --- --- --- %%% %%% --- --- --- %%%  Feedback specific steps:
%%% --- --- --- %%% %%% --- --- --- %%% %%% --- --- --- %%% %%% --- --- --- %%% %%% --- --- --- %%% 
% These stages are skipped if they are turned off. 
    
%%
%%% --- --- --- %%% Just hold in ODT:
    % This will give the trap lifetime, and heating effects if the trap powers are kept at the evaporation end points. 
if Hold_In_ODT == 1
    sq.delay(opt.params);
end 

%% 
%%% --- --- --- %%% Adiabatically ramp up the laser powers:
    %%% Notes:
        % To trap the atoms, we need to ramp the powers up. To determine the ramping 
        % time look at 'makeEquilibriumTime.m', and specify the number of atoms, powers 
        % pre-ramp, and temperature. 

if Ramp_Up_Powers == 1
    %%% Define the number of atoms, and temp:
    N_cond = 2e5;
    T_cond = 100e-9;
    
    %%% Calculate the ramping time:
    TRamp = 2 * get_equilibrium_time(N_cond, T_cond, 2 * opt.keopsys, opt.redpower);
    
    %%% Define a time vector for optical evap: 
    t = linspace(0,TRamp,20);
    
    %%% Define the ODT params from 'opt' 
    final_dipole_ramp.FA = opt.param1; % 1.0 (when halved)
    final_dipole_ramp.RP = opt.param2; % 2.0
        
    %%% Set ODT params:
    sq.find('RedPower CW').after(t,sq.linramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole_ramp.RP))); 
    sq.find('Keopsys FA' ).after(t,sq.linramp(t,sq.find('Keopsys FA' ).values(end),DipolePtoV('keopsys', final_dipole_ramp.FA))); 
    
    %%% Delay for duration of sequence:
    sq.delay(TRamp);
    % sq.delay(opt.params);
end 

% Add section where the power is modulated... 
 
 
%% 
%%% --- --- --- %%% Shadowgraph imaging:
    %%% Notes: 
        % This section takes the images for shadowgraph imaging. It follows on from 
        % the previous section where the dummy images are taken. 
if Shadowgraph_Imaging == 1
    %%% Set the NDI params from 'opt':
    makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_delay,'cycle time',opt.nd.cycle_time,...
         'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_power,'species',85,'num_images',opt.nd.num_images,...
         'pulse delay',opt.nd.pulse_delay);
end

%% 
%%% --- --- --- %%% Thermalisation time:
    %%% Notes: 
        % After taking the images, the cloud is out of equilibrium. Wait another time 
        % determined by 'makeEquilibriumTime.m', now for the post-ramp powers. 

if Thermalisation_Time == 1
    TThermalisation = 250e-3;
    
    %%% Delay for duration of sequence:
    sq.delay(TThermalisation);
end 
 
%% 
%%% --- --- --- %%% Drop those atoms:
timeAtDrop = sq.time;
sq.find('Earth Bias 1').set(0);
sq.find('Earth Bias 2').set(0);
sq.find('Earth Bias 3').set(0);
sq.find('2D MOT Coils').set(0);
sq.find('3DMOT').set(0);
sq.find('87 repump').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(dBtoV('normal',0)); 
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

%%
%%% --- --- --- %%% Information Displayed
if  Display_info == 1
    if CMOT == 1
        if PGC == 1
            if Load_Mag_Evap == 1
                if RF_Evaporation == 1
                    if Optical_Evaporation == 1
                        cprintf('Keywords','Optical evaporation Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                        opt.stage = 'Optical evaporation';
                        baseOpt = evalin('base', 'opt');
                        baseOpt.stage = opt.stage;
                        assignin('base', 'opt', baseOpt);
                    else
                        cprintf('Keywords','RF evaporation Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                        opt.stage = 'RF evaporation';
                        baseOpt = evalin('base', 'opt');
                        baseOpt.stage = opt.stage;
                        assignin('base', 'opt', baseOpt);
                    end
                else
                    cprintf('Keywords','Load into Mag Trap Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                    opt.stage = 'Load into Mag trap';
                    baseOpt = evalin('base', 'opt');
                    baseOpt.stage = opt.stage;
                    assignin('base', 'opt', baseOpt);
                end
            else
                cprintf('Keywords','PGC Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
                opt.stage = 'PGC';
                baseOpt = evalin('base', 'opt');
                baseOpt.stage = opt.stage;
                assignin('base', 'opt', baseOpt);
            end
        else
            cprintf('Keywords','CMOT Stage | Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
            opt.stage = 'CMOT';
            baseOpt = evalin('base', 'opt');
            baseOpt.stage = opt.stage;
            assignin('base', 'opt', baseOpt);
        end
    else
        cprintf('Keywords','MOT Stage |Drop time: %.2f ms | Detuning: %.1f MHz | Camera: %s\n', opt.tof*1000,opt.detuning,Camera)
        opt.stage = 'MOT';
    end

%%
%%% --- --- --- %%% Absorption imaging sequence
    sq.anchor(timeAtDrop);
    if opt.nd.ref_images == 0
        sq.camDelay = timeAtDrop - 2;
    end

    %%% Set the absorption imaging sequence:
    makeImagingSequence(sq,'tof',opt.tof,'pulse time',40e-6,'repump delay',100e-6,...
        'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
        'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
        'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal'); %pulse time use to be 100e-6
    
    %%% turn off anything else hat is left on. 
    sq.find('Redpower CW').set(0);
    sq.find('Redpower TTL').after(100e-6,0);
    sq.find('Keopsys FA').set(0);
    sq.find('Keopsys MO').after(100e-6,0);
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
    sq.find('3DMOT').set(0);
    sq.find('87 repump amp').set(0);
    sq.find('CD0 Fast').set(0);
    sq.waitFromLatest(0.25);
    setSafeValues(sq);
    %%% Compile the results, and run the script. 
    if nargout == 0
        r = RemoteControl;
        r.upload(sq.compile);
        r.run;
    else
        varargout{1} = sq;
    end

end %function