function varargout = makeSequence_Feedback_Jordan_Working(varargin)   

%%% Parse input arguments
    opt = SequenceOptions('load_time',15,'detuning',0,'tof',20e-3,'redpower',2,...
        'raycus',2);
    
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

%%% Initialise the imaging light:
    ImageFreq = (opt.detuning+43.5)^4*6e-8-(opt.detuning+43.5)^3*8e-7-(opt.detuning+43.5)^2*0.0004+(opt.detuning+43.5)*0.0796+5.4703; %assuming imaging field is set to 5V
    dipole_field = 5;
    imaging_field = 0.5;
    ImageAmp = 7; % 9

%%% Initialize sequence
    sq = initSequence;  % load default values (OLD MOT values are default) 
    sq.find('87 imag freq').set(8.3); % RT(ImageFreq)
    sq.find('87 imag amp' ).set(8);

% ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ ~~~ %    

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% The MOT load, and imaging is always on:
Always_on = 1;

%%% Imaging options:
% camera_type = 'low res';
camera_type = 'high res';
% camera_type = 'NDI';

%%% General cooling sequence
CMOT     = 1;
PGC      = 1;
Opt_Pump = 1;
Mag_Trap = 1;
RF_Knife = 1;
ODT_Load = 1;
ODT_Evap = 1;

%%% Transfer to raycus-raycus trap or ramp up powers:
ODT_ramp     = 0; 
ODT_transfer = 1;

%%% Other sequence options:
NDI      = 1;
Therm    = 0;

%%% Type of feedback (if 0 == default)
FB_AOM = 1;

%%% Type of perturbatiopn:
PowerModulate_ODT = 1;
DriveDynamics_MOT = 0;

%%% Include for analysis:
RandomPhase       = 0; % replaces the thermalisation period (turn it off)

%%% Holding in traps for lifetimes:
MT_hold           = 0;
ODT_hold_postRamp = 0;

%%% Optimise:
optimiseAOM_VCO    = 0; % OPTIMISES VCO position. This kicks the cloud at 40 ms into imaging to optimise the VCO bias. Set P_FBL = 0.5 for all stages before the optimise o/w P_FBL = 1.0. 
optimiseAOM_VCO_d  = 0; % OPTIMISES VCO position. Here we drop the atoms for 500 us TOF.

%%% Other things:
TurnOnDipoleBeams = 1; % This should be on if cooling if ODT stages on | turns on during mag trap

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Load into the MOT:
if Always_on == 1

    %%% MOT loading
    sq.find('2DMOT').set(1);
    sq.find('3DMOT').set(1);
    sq.find('87 push').set(1);
    
    %%% 3D MOT beam settings
    sq.find('3DMOT Freq').set(FtoV('trap',26.0)); % 26
    sq.find('3DMOT amp').set(TrapPtoV('trap',1)); % WAS 0.4
    
    %%% 3D repump beam settings
    sq.find('87 repump').set(1);
    sq.find('Repump Switch').set(0);
    sq.find('87 repump freq').set(FtoV('repump',0)); % 0
    sq.find('87 repump amp').set(TrapPtoV('repump',1)); % 0.993
    
    %%% 3D coil settings
    sq.find('H-Bridge Quad').set(1);
    sq.find('CD bit 0').set(0); 
    sq.find('CD bit 1').set(1);
    sq.find('CD2').set(dBtoV('normal',24)); % 25
    sq.find('CD Fine/Fast').set(dBtoV('fine',0)); 
    
    %%% Bias coil settings
    sq.find('Earth Bias 1').set(3); 
    sq.find('Earth Bias 2').set(5); 
    sq.find('Earth Bias 3').set(0); 
    
    %%% Wait for MOT sequence to finish 
    sq.delay(opt.load_time);
        
    %%% Turn off the 2D MOT and coils as well as the push beam
    sq.find('2D MOT Coils').before(10e-3,0);
    sq.find('2DMOT').before(10e-3,0);
    sq.find('87 push').before(10e-3,0);
    sq.find('85 push').before(10e-3,0);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% CMOT sequence
if CMOT == 1

    %%% Time vector:
    Tcmot = 30e-3;
    t = 0 : 2e-3 : Tcmot;
    
    %%% 3D Coils
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(dBtoV('normal',0));
    sq.find('CD Fine/Fast').set(dBtoV('fine',10)); 
    
    %%% Trapping light
    sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',53.4)));
    sq.find('3DMOT amp').set(TrapPtoV('trap',1.45)); % WAS 1.0
    
    %%% Repump
    sq.find('87 repump freq').set(FtoV('repump',12)); % -10.6
    sq.find('87 repump amp').set(TrapPtoV('repump',0.05)); % 0.92
    
    %%% Wait for CMOT sequence to finish
    sq.delay(Tcmot);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% PGC sequence
if PGC == 1
    
    %%% Time vector:
    Tpgc = 6e-3; % 5e-3
    t = linspace(0,Tpgc,10);

    %%% Bias fields:
    sq.find('Earth Bias 1').set(1.2);  % Probably N/S % Originally OFF
    sq.find('Earth Bias 2').set(6.5);  % Probably U/D % Originally OFF
    sq.find('Earth Bias 3').set(0.05); % Probably E/W % Originally OFF

    %%% Turn off magnetic fields:
    sq.find('CD fine/fast').set(dBtoV('fine',0));
    sq.find('CD0 Fast').set(dBtoV('normal',0));
    
    %%% Change the trap frequency and amplitude
    sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',79.8)));  % RT(76.8)
    sq.find('3DMOT amp' ).after(t,sq.minjerk(t,sq.find('3DMOT amp' ).values(end),TrapPtoV('trap',0.85))); % WAS 0.5
    
    %%% Change the rempump frequency and amplitude
    sq.find('87 repump freq').set(FtoV('repump',11)); % 11
    sq.find('87 repump amp').set(TrapPtoV('repump',0.05)); %0.25
    
    %%% Wait for PGC to finish
    sq.delay(Tpgc);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Optical pump atoms into the F = 1 manifold
if Opt_Pump == 1
    
    %%% Time vector:
    Tdepump = 1e-3;

    %%% Turn the repump light off:
    sq.find('repump switch').set(1); %fiber switch off (it's inverted)
    sq.find('87 repump').set(0);
    sq.find('87 repump amp').set(0);
    
    %%% Wait for pumping to finish:
    sq.delay(Tdepump);

    %%% Turn the MOT off:
    sq.find('3DMOT').set(0);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Load into magnetic trap
if Mag_Trap == 1

    %%% %%% %%% Magnetic loading/trapping:
        %%% Time vector
        Tmagload = 150e-3; % 200e-3
        t = 0 : 10e-3 : Tmagload;
    
        %%% Turn the earth cancelling biases off:
        sq.find('Earth Bias 1').set(0); % Probably N/S
        sq.find('Earth Bias 2').set(0); % Probably U/D
        sq.find('Earth Bias 3').set(0); % Probably E/W
    
        %%% Turn on the magnetic trap:
        dBLoad = 110;
        sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
        sq.find('CD Fine/Fast').set(dBtoV('fine', 0)); % 4
        % sq.delay(Tmagload); % just holding in MT
    
    %%% %%% %%% Optical loading/trapping:
        if TurnOnDipoleBeams == 1        
            %%% Load into optical dipole trap
            Toptload = 400e-3; % 100e-3
            t = 0 : 20e-3 : Toptload;

            %%% Turn on the optical dipole traps
            sq.find('Raycus TTL').set(1); % 3.9
            sq.find('Raycus CW' ).after(t,sq.minjerk(t,0, DipolePtoV('Raycus', 8))); % WAS 13
            sq.find('Redpower TTL').set(1); % RT(1)
            sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower', 13))); % WAS 15
            
            %%% Turn on FBL if using AOM deflection
            if FB_AOM == 1
                sq.find('feedback laser ttl').set(1);
                sq.find('feedback laser power').after(t,sq.minjerk(t,0,DipolePtoV('feedback', 0.5))); % WAS 0.5
            end 
        end 
        
        %%% Turn off the MOT beams:
        sq.find('MOT bias').set(1);
        sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));

        %%% Wait for magnetic trapping and optical trapping to finish:
        if TurnOnDipoleBeams == 1
            sq.delay(max(Toptload,Tmagload));
        else 
            sq.delay(Tmagload);
        end 
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Hold in magnetic trap or dimple trap
if MT_hold == 1
    Thold = 1; % WHEN CHARACTERISING
    sq.delay(Thold);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% RF evaporation
if RF_Knife == 1
    
    %%% Time vector from evap params:
    rf_start = 16.0; % 16.0
    rf_end   = 2.0;  % 2.0
    rf_rate  = 2.0;  % 2.0
    Tevap    = (rf_start - rf_end)/rf_rate;
    rf_ramp_type = 'lin';
    rf_exp_time_constant = 2;
    t = linspace(0,Tevap,51);
    
    %%% Ramp the RF frequency down:
    sq.find('RF atten').set(1);
    if strcmpi(rf_ramp_type,'exp')
        sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
    elseif strcmpi(rf_ramp_type,'lin')
        sq.find('RF frequency').set(RFtoV(rf_start));
        sq.delay(0.25);
        sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
    end

    %%% Delay for RF evap to finish:
    sq.delay(Tevap);
    
    %%% Turn RF field off:
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Take dummy images
if NDI == 1
    time_at_evap_end = sq.time;
    %%% If we are taking images we take images -- shocking right
    if opt.nd.ref_images > 0       
        sq.anchor(sq.time - 3);
        sq.camDelay = sq.time - 2;
        
        %%% Turn on FBL for image subtraction
        if FB_AOM == 1
            P_FBL = 1.2;
            sq.find('Feedback laser power').set(DipolePtoV('Feedback', P_FBL)); % Set FB laser to 1.2 W 
        end 

        makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',5e-6,'cycle time',100e-3,...
            'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_power,'species',85,'num_images',opt.nd.ref_images,...
            'pulse delay',opt.nd.pulse_delay);      
    end
    
%     %%% Vary the 85 image freq and amplitude for PA work:
%     if opt.nd.ref_images > 0
%         time_at_evap_end = sq.time;
%         sq.anchor(sq.time - 3);
%         sq.camDelay = sq.time - 2;
% 
%         %%% Set the frequency and amp to get 10 mW:
%         v85Freq = opt.param1;
%         v85Amp  = opt.param2;
%         v85Amp  = TrapPtoV_NDI(0.1e-3); % low power to avoid heating
% 
%         makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',5e-6,'cycle time',100e-3,...
%             'imaging freq',v85Freq,'imaging voltage',v85Amp,'species',85,'num_images',opt.nd.ref_images,...
%             'pulse delay',opt.nd.pulse_delay);
%         sq.anchor(time_at_evap_end);
%      end

    sq.anchor(time_at_evap_end);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Turn off magnetic trap
if ODT_Load == 1
%       %%% OLD THINGS: (RT does not use this)   
%     Trampcoils = 0.4; % RT(0.8)
%     dB_weak = 0;
%     t = linspace(0,Trampcoils,51);
%     sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));
%     sq.delay(Trampcoils);

    %%% Set the dipole beam back to nominal value:
    if FB_AOM == 1
        P_FBL = 1.00;
        sq.find('Feedback laser power').set(DipolePtoV('Feedback', P_FBL)); % Set FB laser to 1.0 W 
    end 

    %%% Turn off the magnetic trap
    sq.find('CD0 Fast').set(0);
    sq.find('CD Fine/Fast').set(0);

    %%% Hold times:
    if(0)
        %%% Turn of lasers:
        % sq.find('Raycus CW'   ).set(0);
        % sq.find('Raycus TTL'  ).set(0);
        % sq.find('Redpower CW' ).set(0);
        % sq.find('Redpower TTL').set(0);

        %%% Add delay:
        % sq.delay(50e-3); % IF HOLDING IN SINLGE DIPOLE BEAMS USE 50e-3;
        % sq.delay(1); % HOLD for 1s when benchmarking
    end 
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Optical evaporation
if ODT_Evap == 1
    
    %%% Time vector:
    Tevap = 3.5; % 3.0
    Tc = 0.5;
    t = linspace(0,Tevap,51);

    %%% Set the evaporation end points of the trapping beams:
    final_dipole.RP = opt.redpower; %  || 1.3 (NEW)
    final_dipole.FA = opt.raycus;   %  || 1.2 (NEW)
    sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),Tc)); % 0.50
    sq.find('Raycus CW'  ).after(t,sq.expramp(t,sq.find('Raycus CW'  ).values(end),DipolePtoV('raycus',  final_dipole.FA),Tc)); % 0.50

    %%% CHECKING OVERLAP
    if optimiseAOM_VCO_d == 1
        final_dipole.FB = 1.5;   
        sq.find('feedback laser power').after(t,sq.expramp(t,sq.find('feedback laser power').values(end),DipolePtoV('feedback', final_dipole.FB),0.5)); % 0.50
    end 
    
    %%% Wait for the optical evaporation to finish:
    sq.delay(Tevap);

    %%% Delays for checking overlaps:
    % sq.find('Raycus CW'  ).after(t,sq.expramp(t,sq.find('Raycus CW'  ).values(end),DipolePtoV('raycus',  0),0.3)); % 
    % sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',0),0.3)); % 
    % sq.delay(0.2); % x-position
    % sq.delay(1.0); % z-position
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Hold in ODT (pre ramp):
if optimiseAOM_VCO_d == 1
    Tramp = 0.50;
    Thold = 0.25; % FOR CHARACTERISATION
    
    %%% Check if beams are overlapped
    check_overlap.RP = 0;  
    check_overlap.FA = 0;  
    check_overlap.FB = 1.5;  
    sq.find('RedPower CW'         ).after(t,sq.expramp(t,sq.find('RedPower CW'         ).values(end),DipolePtoV('redpower', check_overlap.RP),Tramp)); % 0.50
    sq.find('Raycus CW'           ).after(t,sq.expramp(t,sq.find('Raycus CW'           ).values(end),DipolePtoV('raycus'  , check_overlap.FA),Tramp)); % 0.50
    sq.find('feedback laser power').after(t,sq.expramp(t,sq.find('feedback laser power').values(end),DipolePtoV('feedback', check_overlap.FB),Tramp)); % 0.50
    
    %%% Add delay
    sq.delay(Thold);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%

% %%% Adiabatically ramp up trap powers:
% if ODT_ramp == 1
%     %%% Set ODT powers to ramp to
%     final_dipole_ramp.FA  = opt.raycus   + 0.1; % NORMAL USE
%     final_dipole_ramp.RP  = opt.redpower + 0.0; % NORMAL USE 
%     if optimiseAOM_VCO == 1
%         final_dipole_ramp.FBL = 0.5;
%         sq.find('feedback laser power').set(DipolePtoV('feedback',final_dipole_ramp.FBL)); % Lower power before kick
%     else 
%         final_dipole_ramp.FBL = 1.0;
%         sq.find('feedback laser power').set(DipolePtoV('feedback',final_dipole_ramp.FBL)); % NORMAL USE
%     end 
% 
%     %%% Time vector:
%     TAdiabaticRamp_up = 200e-3; % NORMAL 200e-3 
%     t = linspace(0,TAdiabaticRamp_up,51);
%        
%     %%% Set ODT params:
%     sq.find('RedPower CW'         ).after(t,sq.linramp(t,sq.find('RedPower CW'         ).values(end),DipolePtoV('redpower',final_dipole_ramp.RP))); 
%     sq.find('Raycus CW'           ).after(t,sq.linramp(t,sq.find('Raycus CW'           ).values(end),DipolePtoV('raycus',  final_dipole_ramp.FA))); 
%     sq.find('feedback laser power').after(t,sq.linramp(t,sq.find('feedback laser power').values(end),DipolePtoV('feedback',final_dipole_ramp.FBL)));
%     
%     %%% Wait for the ramping to finish:
%     sq.delay(TAdiabaticRamp_up);
% end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Hold in ODT (post ramp):
if ODT_hold_postRamp == 1
    Thold = 1;
    sq.delay(Thold);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Adiabatically transfer atoms:
if ODT_transfer == 1
    %%% Set ODT powers to ramp to
    final_dipole_ramp.FA  = opt.raycus + 0.4; % RAMP up by 0.6
    final_dipole_ramp.RP  = 0.1; % RAMP down
    final_dipole_ramp.FBL = 1.2; 

    %%% Time vector
    TAdiabaticRamp = 100e-3;  
    TDelay         = 0;
   
    t = linspace(0,TAdiabaticRamp,51);

    %%% Apply the ramps:
    sq.find('Raycus CW'           ).after(t,sq.linramp(t,sq.find('Raycus CW'           ).values(end),DipolePtoV('raycus',   final_dipole_ramp.FA ))); 
    sq.find('feedback laser power').after(t,sq.linramp(t,sq.find('feedback laser power').values(end),DipolePtoV('feedback', final_dipole_ramp.FBL))); 

    %%% Add delay between ramping lasers:
    sq.delay(TDelay);

    %%% Ramp down redpower
    sq.find('RedPower CW'         ).after(t,sq.linramp(t,sq.find('RedPower CW'         ).values(end),DipolePtoV('redpower',final_dipole_ramp.RP ))); 
    
    %%% Add delay:
    sq.delay(TAdiabaticRamp);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Power modulation of ODT for breathing mode:
if PowerModulate_ODT == 1
    %%% Define the modulation parameters to time and cycles:
    nn = 5; % 10 if FBL o/w 4
    f_mod = sqrt(5/2) * 21.5; % 21????
    T_mod = 1/f_mod;
    dt = T_mod * 0.05; % discretisation time
    t = 0 : dt : nn * T_mod;

%     %%% Modulate powre of trapping laser to induce mode:
%     if(0)
%         %%% Define the modulation to power:
%         final_dipole_ramp.FA = opt.raycus;
%         final_dipole_ramp.RP = opt.redpower;
%         dP_red = 0.15; % 0.25
%         dP_ray = 0.00; % 0.00
%         final_dipole_mod.RP = final_dipole_ramp.RP * (1 + dP_red * sin(2 * pi * f_mod * t)); 
%         final_dipole_mod.FA = final_dipole_ramp.FA * (1 + dP_ray * sin(2 * pi * f_mod * t));
%         sq.find('RedPower CW').after(t, DipolePtoV('redpower', final_dipole_mod.RP));
%         sq.find('Raycus CW'  ).after(t, DipolePtoV('raycus'  , final_dipole_mod.FA));
% 
%         %%% Add delays:
%         sq.delay(nn * T_mod); % accounts for modulation
%     end 

    %%% Modulate power of FBL to induce mode:
    if(1)
        %%% Define the modulation to power:
        P_FBL = 1.2;
        dP_mod = 0.50; % 50% modulation depth
        final_dipole_mod.FBL = P_FBL * (1.0 + dP_mod * sin(2 * pi * f_mod * t));
        sq.find('feedback laser power').after(t, DipolePtoV('feedback', final_dipole_mod.FBL));             
        
        %%% Add delays:
        sq.delay(nn * T_mod); % accounts for modulation
    end 
end 

%%% Set FBL power to operational power:
% sq.find('feedback laser power').set(DipolePtoV('feedback',final_dipole_ramp.FBL)); % NORMAL USE

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 

% %%% Ramp the MOT coils for vertical motion:
% if DriveDynamics_MOT == 1
%     %%% Define parameters :
%     T_MOT    = 100e-3;  % Time that coils are on
%     dB_MOT   = 5.0;    % Offset from unit step
%     % dB_rand  = 0.0;    % Gets multiplied by dB_drive
%     dB_drive = 8.0;    % Driving amplitude (must less than offset for first driving configuration)
%     
%     %%% Time
%     dt    = T_MOT * 0.01;   % Discretisation time
%     t     = 0 : dt : T_MOT; % Time vector
% 
%     %%% Set MOT coil params:
%         %%% Unit impulse
%         if(0)
%             sq.find('CD Fine/Fast').set(dBtoV('fine', dB_MOT)); % Sets the field to dB_MOT
%         end 
%     
%         %%% Continuous driving
%         if(0)
%             %%% Set a seed (which can be used to construct driving profiles later)
%             % rng(1); %%% Replace with opt.params for different seeds
% 
%             %%% Define the driving field (at resonance with x mode)
%             fx     = sqrt(5/2) * 22.5; % min(sort(get_trap_freq(final_dipole_ramp.FA, final_dipole_ramp.RP)));
%             % dB_MOT = dB_drive * sin(2 * pi * fx / 2 * t).^2 .* (1 + dB_rand * randn(1, length(t)));
%             dB_MOT = (1 - exp(-30  * t)) .* dB_drive .* sin(2 * pi * fx * t / 2) .^2;
% 
%             %%% Set magnetic field gradient
%             sq.find('CD Fine/Fast').after(t, dBtoV('fine', dB_MOT)); % 
%         end 
% 
%     %%% Add delays:
%     sq.delay(T_MOT); % accounts for modulation
%     
%     %%% Set coils back to zero:
%     sq.find('CD Fine/Fast').set(dBtoV('fine', 0)); % 
% end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Just kick the atoms in trap
if optimiseAOM_VCO == 1
    P_FBL = 1.50;
    sq.find('feedback laser power').set(DipolePtoV('feedback',P_FBL)); % NORMAL USE
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Non-destructive imaging
if NDI == 1
    %%% Previously used things by RT:
% num_mod_cycles = 4;
% fmod = 32;
% Tmod = 1/fmod;
% dt = 0.05*Tmod;
% t = 0:dt:(num_mod_cycles*Tmod);
% sq.find('RedPower CW').after(t,DipolePtoV('redpower',1.34 + 0.25*sin(2*pi*fmod*t)));
% sq.find('Raycus CW').after(t,DipolePtoV('raycus',final_dipole.FA - 0.125/2.25*sin(2*pi*fmod*t)));
% sq.delay(num_mod_cycles*Tmod);
% % sq.delay(opt.params(1));
% sq.delay(10e-3);
% sq.find('Raycus CW').set(DipolePtoV('raycus',final_dipole.FA + 0.0));
% sq.find('RedPower CW').after(50e-3,DipolePtoV('redpower',final_dipole.RP + 0.0));
% Tdrop = 0.5e-3;
% sq.find('Raycus CW').set(0).after(Tdrop,sq.find('Raycus CW').values(end-1));
% sq.delay(Tdrop + 0.5e-3);
% Tpulse = 15e-3;
% sq.find('RedPower CW').set(sq.find('RedPower CW').values(end) + 0.2).after(Tpulse,sq.find('RedPower CW').values(end-1));
% sq.delay(Tpulse);
% sq.delay(500e-3);
% sq.anchor(sq.time - num_mod_cycles*Tmod/2);
    
    %%% If you want to kick the cloud to optimise AOM VCO
    if optimiseAOM_VCO == 1
        sq.anchor(sq.time - 40e-3);
    end 

    %%% Turn on the FB laser to near max power:
    if optimiseAOM_VCO == 1
        
    else
        P_FBL = 1.0;
        sq.find('feedback laser power').set(DipolePtoV('feedback',P_FBL)); % 
    end 

    %%% NDI sequence:
    makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_delay,'cycle time',opt.nd.cycle_time,...
        'imaging freq',8.5,'imaging amplitude',opt.nd.pulse_power,'species',85,'num_images',opt.nd.num_images,...
        'pulse delay',opt.nd.pulse_delay);

%     %%% Set the frequency and amp to get 10 mW:
%     v85Freq = opt.param1;
%     v85Amp  = opt.param2;
%     
%     makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_delay,'cycle time',opt.nd.cycle_time,...
%         'imaging freq',v85Freq,'imaging voltage',v85Amp,'species',85,'num_images',opt.nd.num_images,...
%         'pulse delay',opt.nd.pulse_delay); 
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Thermalisation
if Therm == 1
    %%% Wait for sample to thermalise
    TThermalisation = 200e-3; % WAS 250e-3
    sq.delay(TThermalisation);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Random Phase delay
if RandomPhase == 1
    %%% If adding a random amount of time 
    % TRanddomPhase = opt.params;
    TRanddomPhase = randi([150, 250])*1e-3;
    sq.delay(TRanddomPhase);
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%%
%%% Drop atoms
if Always_on == 1
% sq.delay(1);
    timeAtDrop = sq.time;
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
    sq.find('Raycus CW').set(0);
    sq.find('Raycus TTL').after(100e-6,0);
    sq.find('feedback laser power').set(0);
    sq.find('feedback laser ttl').after(100e-6,0);
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
    % sq.find('Variable wave plate').set(-3.7);   %This value switches to absorption imaging
end 

%%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% %%% 
%%% Take Absorption Image
if Always_on == 1
    
    sq.anchor(timeAtDrop);
    if opt.nd.ref_images == 0
        sq.camDelay = timeAtDrop - 2;
    end
    
    if strcmpi(camera_type,'high res')
    %%% Main imaging sequence (regular camera):
    makeImagingSequence(sq,'tof',opt.tof,'pulse time',80e-6,'repump delay',100e-6,...
        'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
        'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
        'fibre switch delay',1e-3,'imaging_field',dipole_field,'image type','horizontal');
    elseif strcmpi(camera_type,'low res')
    %%% Low resolution camera:
    makeImagingSequence(sq,'tof',opt.tof,'pulse time',40e-6,'repump delay',100e-6,...
        'repump time',200e-6,'cam time',5e-6,'cycle time',500e-3,...
        'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
        'fibre switch delay',1e-3,'imaging_field',dipole_field,'image type','horizontal');
    elseif strcmpi(camera_type,'NDI')
    makeImagingSequence(sq,'tof',opt.tof,'pulse time',80e-6,'repump delay',100e-6,...
        'repump time',200e-6,'cam time',5e-6,'cycle time',200e-3,...
        'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
        'fibre switch delay',1e-3,'imaging_field',dipole_field,'image type','85');
    end 
    
    % turn off the dipoles
    sq.find('Redpower CW').set(0);
    sq.find('Redpower TTL').after(100e-6,0);
    sq.find('Raycus CW').set(0);
    sq.find('Raycus TTL').after(100e-6,0);
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(RFtoV(20));
    sq.find('3DMOT').set(0);
    sq.find('87 repump amp').set(0);

    sq.waitFromLatest(0.25);
    setSafeValues(sq);
 
    if nargout == 0
        r = RemoteControl;
        r.upload(sq.compile);
        r.run;
    else
        varargout{1} = sq;
    end
end 

%%% %%% %%% %%% %%% Ends script
end