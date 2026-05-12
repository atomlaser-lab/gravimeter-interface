function varargout = makeSequenceGinger(varargin)
%% Parse input arguments
opt = parse_maker_variable_argument_list(varargin{:});

% %make SLM pattern
l = opt.param3; %LG charge
% dx = round(800.5*0 - 0.008005);
dx = 4;
% dy = round(400.3*0 - 0.004003);
dy = -1;
LG_SPP_generator(l,dx,dy);
% LG_SPP_generator_rot(l, dx, dy, 0, false);
% LG_SPP_generator_general(l, dx, dy, 1.5, 0)

% % % % for HG
% n = 1;
% m = 0;
% HG_SPP_generator(m, n, dx, dy, 256, 22)

% % for saddle phase
% Saddle_SPP_generator(l, dx, dy)

ImageFreq = opt.detuning + 0.5; %Low intensity after dipole evaporation, high res imaging
% ImageFreq = opt.detuning + 1; %Low intensity after dipole evaporation, vertical imaging
dipole_field = 1; %In Gauss
ImageAmp = 0.1;
%% Initialize sequence
sq = initSequence;  %load default values (OLD MOT values are default)
sq.find('87 imag freq').set(ImageFreq);
sq.find('87 imag amp').set(1);

if opt.stage.use_dipoles
    sq.find('Raycus TTL').set(1);
    sq.find('Raycus CW').set(200e-3);
    sq.find('RedPower TTL').set(1);
    sq.find('RedPower CW').set(100e-3);
end

%% MOT loading
if opt.stage.use_mot
    sq.find('Vertical MOT Mirror').set(0);
    sq.find("3DMOT").set(0);
    sq.delay(0.5);
    sq.find('2DMOT Freq').set(18);
    sq.find('Push Freq').set(5);
    %     sq.find('Push amp').set(3.75);
    sq.find('2DMOT').set(1);
    sq.find('3DMOT').set(1);
    sq.find('87 push').set(1);
    % 3D MOT beam settings
    sq.find('3DMOT Freq').set(18);
    sq.find('3DMOT amp').set(1);
    % 3D repump beam settings
    sq.find('87 repump').set(1);
    sq.find('Repump shutter').set(1);
    %     sq.find('Repump Switch').set(0);
    sq.find('87 repump freq').set(0);
    sq.find('87 repump amp').set(1);
    % 3D coil settings
    sq.find('H-Bridge Helm').set(0); %Quad is the default setting
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(14); %Coarse control of 3D coils %14
    sq.find('CD Fine/Fast').set(0); % fine control of 3D coils
    % Bias coil settings
    sq.find('Bias E/W').set(0.4);
    sq.find('Bias N/S').set(3);
    sq.find('Bias U/D').set(6);
    %Delay for the load time
    sq.delay(opt.load_time);
    %
    % Turn off the 2D MOT and coils as well as the push beam
    %
    sq.find('2D MOT Coils').before(10e-3,0); %active low
    sq.find('2DMOT').before(10e-3,0);
    sq.find('87 push').before(10e-3,0);
    sq.find('C6 - N/C').set(1).after(100e-6,0); %trigger SLM
end
%% CMOT sequence
%
% Apply a compressed MOT sequence to temporarily increase the density by
% reducing spontaneous emission.  We switch to CD channel 0b00 = 0 because
% it is the fast channel
%
if opt.stage.use_cmot
    Tcmot = 5e-3;
    t = 0:0.25e-3:Tcmot;
    % 3D Coils
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(0);
    sq.find('CD Fine/Fast').set(0.5);
    %Trapping light
    sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),55));
    sq.find('3DMOT amp').set(1);
    %Repump
    sq.find('87 repump freq').set(2.5); %-7
    sq.find('87 repump amp').set(1);

    sq.delay(Tcmot);
end
%% PGC sequence
%
% Apply polarization gradient cooling to reduce the temperature of the
% atoms.  We use CD channel 0b00 = 0 as it is the fast-switching channel
%
if opt.stage.use_pgc
    Tpgc = 2e-3;
    t = 0:0.25e-3:Tpgc;
    % t = linspace(0,Tpgc,26);
    %     sq.find('Bias E/W').set(0);
    %     sq.find('Bias N/S').set(0);
    sq.find('CD fine/fast').set(0);
    sq.find('CD0 Fast').set(0);
    sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),75));
    sq.find('3DMOT amp').after(t,sq.linramp(t,sq.find('3DMOT amp').values(end),0.9)); %0.5

    sq.find('87 repump freq').set(-5);%-4.8
    sq.find('87 repump amp').set(1);%0.004

    sq.delay(Tpgc);
end

%% Optical pump atoms into the F = 1 manifold
%
% Turn off repump field so that atoms are optically pumped into the F = 1
% manifold.
%
if opt.stage.use_pump
    Tdepump = 3e-3;
    sq.find('Repump shutter').set(0);
    sq.find('87 repump').set(0).after(5e-3,1);
    %     sq.find('87 repump').set(0);
    sq.find('87 repump freq').set(20);
    sq.find('3DMOT freq').set(75);
    sq.delay(Tdepump);
    sq.find('3DMOT').set(0);
end

%% Load into magnetic trap
%
% Load into the magnetic trap at a high gradient.  We switch quickly to a
% low value and then ramp up to the target value
%
if opt.stage.use_mag
    Tmagload = 150e-3;
    t = 0:10e-3:Tmagload;
    dBmax = 110; %110
    dBLoad = 55;
    sq.find('CD0 Fast').after(t,sq.linramp(t,dBLoad,dBmax));
    sq.find('CD Fine/Fast').set(0);
    sq.find('Vertical MOT Mirror').set(1);
    sq.delay(Tmagload);

    if opt.stage.use_dipoles
        Toptload = 400e-3;
        t = 0:20e-3:Toptload;
        sq.find('Raycus TTL').set(1);
        sq.find('Redpower TTL').set(1);
        sq.find('Raycus CW').after(t,sq.linramp(t,0,4));
        sq.find('Redpower CW').after(t,sq.linramp(t,0,12));
        sq.delay(max(Toptload - Tmagload,0));
    end

    if ~opt.stage.evap_mag
        sq.delay(1);
    end
end

%% RF evaporation
%
% Remove hot atoms from the sample using RF transitions between the trapped
% |F = 1, m_F = -1> state and the untrapped |F = 1, m_F = 0> state.  All
% frequencies are in MHz
%
if opt.stage.use_evap_mag
    rf_start = 20; %16
    rf_end =   0.8; %0.8 %0.75
    %         rf_end = 4;
    rf_rate = 3;    %MHz/s 3
    Tevap = (rf_start - rf_end)/rf_rate;
    t = linspace(0,Tevap,50);

    sq.find('DDS switch').set(0);
    sq.find('RF switch').set(1);
    sq.find('RF frequency').set(rf_start);
    sq.delay(0.25);
    sq.find('RF frequency').after(t,sq.linramp(t,rf_start,rf_end));
    sq.delay(Tevap);

    %     sq.find('Repump shutter').set(1); %opening this now to pump into F = 2 in ODT

    sq.find('RF switch').set(0);
    sq.find('RF Frequency').set(20);
end

%% Load into dipole trap
if opt.stage.use_dipoles
    Trampcoils = 0.3; %0.3
    dB_weak = 0;
    t = linspace(0,Trampcoils,51);
    sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dB_weak));
    sq.find('MOT bias coil').after(t,sq.linramp(t,sq.find('MOT bias coil').values(end),dipole_field));
    sq.delay(Trampcoils);
%     sq.delay(1);
end

%% Optical evaporation
if opt.stage.use_evap_dipoles

    %     %pump into F = 2
    %     sq.find('87 repump').set(1);
    %     sq.find('87 repump freq').set(0);
    %     sq.find('87 repump amp').set(1);
    %     sq.delay(5e-3);
    %     sq.find('87 repump amp').set(0);

    %blow away F=2 atoms
    Tblow = 10e-6;
    sq.find('3DMOT').set(1);
    sq.find('3DMOT Freq').set(18); %FtoV('trap',0)
    sq.find('3DMOT amp').set(1); % WAS: 0.4 NOW 0.7
    sq.delay(Tblow);
    sq.find('3DMOT').set(0);

    Tevap = 5;
    TC = 0.75;
    t = linspace(0,Tevap,51);
    sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),opt.raycus,TC)); %opt.redpower
    sq.find('Raycus CW').after(t,sq.expramp(t,sq.find('Raycus CW').values(end),opt.raycus,TC));
    sq.delay(Tevap);
    sq.delay(1);
end

% Trigger the DDS
    sq.ddsTrigDelay = sq.time;
    sq.find('DDS Trigger').before(10e-3,1).after(10e-3,0);%.after(1e-3,1);
%     sq.find('H-Bridge Helm').set(1);

%%% ARP with DDS %%%
if (1)
    sq.find('H-Bridge Helm').set(1);
    sq.delay(0.1);
    sq.find('DDS switch').set(1);

%     % Trigger the DDS
%     %
%     sq.ddsTrigDelay = sq.time;
%     sq.find('DDS Trigger').before(10e-3,1).after(10e-3,0);
%     %
%     % Switch to the Helmholtz configuration
    %
    sq.find('CD0 Fast').set(0);
    sq.delay(10e-3);
    T_helmholtz_ramp = 100e-3;
    t = linspace(0,T_helmholtz_ramp,51);
    sq.find('CD0 Fast').after(t,sq.linramp(t,0,20));
    sq.delay(T_helmholtz_ramp);
    %     sq.delay(100e-3);
    %
    % Apply RF
    %
    sq.dds(1).set(38,0,0); %38.5
    sq.dds(2).set(110,0,0);
    sq.delay(10e-6);
    %     Tarp = 20e-3; %10e-3
    span = 1;
    Tarp = 20e-3/0.5*span;
    t = linspace(0,Tarp,501);
    df = 38 + span*sq.linramp(t,-0.5,0.5);
    w = Tarp/3;
    amp = 0.5*exp(-2*(t - Tarp/2).^2/w^2);
    sq.dds(1).after(t,df,amp,0);
    sq.dds(2).after(t,110,0,0);
    sq.delay(Tarp);
    sq.dds(1).set(110,0,0);
    sq.dds(2).set(110,0,0);
    sq.find('CD0 Fast').set(0.8); %set to 1 for Raman
    sq.find('RF Switch').set(0);
    sq.find('DDS Switch').set(0);
    %     sq.delay(0.1);
    %     sq.find('H-Bridge Helm').set(1);
%     sq.find('Repump shutter').set(1);
    sq.delay(20e-3);
end

%% Drop atoms
timeAtDrop = sq.time;
sq.find('2D MOT Coils').set(0);
sq.find('3DMOT').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(0.8); %set to 1 for Raman
sq.find('MOT bias coil').set(0);
sq.find('CD2').set(0);
sq.find('CD Fine/Fast').set(0);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('DDS Switch').set(0);
sq.find('RF switch').set(0);
sq.find('RF Frequency').set(20);
sq.find('Raycus CW').set(0);
sq.find('Raycus TTL').set(0);
sq.find('RedPower CW').set(0);
sq.find('RedPower TTL').set(0);
sq.find('MOT bias').set(0);
sq.find('Repump shutter').set(0);

% sq.find('Bias E/W').set(0); %0.4
% sq.find('Bias N/S').set(0); %3
% sq.find('Bias U/D').set(0); %6

%% VMG! %%
if(1)
    sq.anchor(timeAtDrop);
    raman_delay_drop = 10e-3; %10e-3
    interrogation_time = opt.param2; %10e-6
    Traman_pi1 = 16e-6/2; %12e-6/2
    Traman_pi2 = Traman_pi1;
%     Traman_pi2 = round((interrogation_time*1e6*0.000572+11.995)/2)*1e-6;     
    F2_tof = raman_delay_drop + interrogation_time + 0.5e-3;
%     F2_tof = raman_delay_drop + 5e-3;
    F1_tof = F2_tof + 3e-3; %3e-3 minimum with smallest frame for BEC
    F1_minus_F2_tof = F1_tof - F2_tof; %needs to be 12ms minimum for full frame
    AI_phi_LG = 0;
    AI_phi_G = opt.param1;
    sq.anchor(timeAtDrop + raman_delay_drop);
    Power_raman_G = 1; %1 %sideband
    Power_raman_LG = 1; %0.1 %carrier
    Delta_raman1 = 19.9764; %19.9764
    Delta_raman2 = Delta_raman1;
%     Delta_raman2 = 0.00000125*interrogation_time*1e6+19.976;

%      if opt.params == 1

%     %     single Raman pulse
%     sq.dds(1).set(110+Delta_raman1/4,Power_raman_LG,0);
%     sq.dds(2).set(110-Delta_raman1/4,Power_raman_G,0);
%     sq.delay(Traman_pi1);
%     sq.dds(1).set(110,0,0);
%     sq.dds(2).set(110,0,0);
%     sq.delay(10e-6);

%      else
    
        %     BS1
        sq.dds(1).set(110+Delta_raman1/4,Power_raman_LG,0); %sideband
        sq.dds(2).set(110-Delta_raman1/4,Power_raman_G,0); %carrier
        sq.delay(Traman_pi1);
        sq.dds(1).set(110+Delta_raman1/4,0,AI_phi_LG);
        sq.dds(2).set(110-Delta_raman1/4,0,AI_phi_G);
        %     sq.find('LG Shutter').set(1).after(100e-6,0); %trigger SLM
    
        sq.delay(interrogation_time);
    
        %     BS2
        sq.dds(1).set(110+Delta_raman2/4,Power_raman_LG,AI_phi_LG); %sideband
        sq.dds(2).set(110-Delta_raman2/4,Power_raman_G,AI_phi_G); %carrier
        sq.delay(Traman_pi2);
        sq.dds(1).set(110+Delta_raman2/4,0,0);
        sq.dds(2).set(110-Delta_raman2/4,0,0);
        sq.delay(10e-6);
%     
%      end

%     %     %turning field off and imaging field on
    sq.find('CD0 Fast').set(0); %0.133 ~0.5G

    trap_hold = 0;
    sq.delay(trap_hold);
    % sq.find('87 imag').set(0);
% 
%     sq.find('Raycus CW').set(0);
% sq.find('Raycus TTL').set(0);
% sq.find('RedPower CW').set(0);
% sq.find('RedPower TTL').set(0);

end

%% Stern-Gerlach

if (0)
    % HEY RYAN! MAKE SURE THIS IS UNCOMMENTED TO SWITCH BACK TO ANTI-HELMHOLTZ CONFIGURATION FOR STERN-GERLACH
    sq.delay(10e-6);
    sq.find('CD0 Fast').set(0);

    sq.anchor(timeAtDrop);
    sq.delay(1e-3); %1e-3
    %     sq.delay(raman_delay_drop+Traman_pi+interrogation_time + 0.5e-3); %1e-3
    SG_Pulse_time = 5e-3;
    SG_Amp = 20; %80
    t = linspace(0,SG_Pulse_time,50);
    sq.find('CD0 Fast').after(t,sq.linramp(t,0,SG_Amp));
    sq.delay(SG_Pulse_time);
    sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),0));
end

%% Take Absorption Image
sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 3;
% sq.delay(trap_hold);
% % sq.delay(raman_delay_drop + Traman_pi + interrogation_time);

% %for horizontal
% % high res
% makeImagingSequence(sq,'tof',opt.tof,'pulse time',4*40e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
%     'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',1,...
%     'repump shutter delay',2e-3,'imaging_field',dipole_field,'image type','horizontal');

makeImagingSequence_4_Images(sq,'tof',F2_tof,'tof2',F1_minus_F2_tof,'pulse time',4*40e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',5e-6,'cycle time',300e-3,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp*10,...
    'repump shutter delay',2e-3,'imaging_field',dipole_field,'image type','horizontal');

% %low res
% makeImagingSequence(sq,'tof',opt.tof,'pulse time',4*40e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
%     'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',10*ImageAmp,...
%     'repump shutter delay',2e-3,'imaging_field',dipole_field,'image type','horizontal');

% makeImagingSequence_4_Images(sq,'tof',F2_tof,'tof2',F1_minus_F2_tof,'pulse time',4*40e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',300e-3,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp*10,...
%     'repump shutter delay',2e-3,'imaging_field',dipole_field,'image type','horizontal');

% % % %for vertical
% makeImagingSequence(sq,'tof',opt.tof,'pulse time',4*40e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',300e-3,...
%     'manifold',2,'imaging freq',ImageFreq,'imaging amplitude',1,...
%     'repump shutter delay',2e-3,'imaging_field',dipole_field,'image type','vertical');

% makeImagingSequence_4_Images(sq,'tof',F2_tof,'tof2',F1_minus_F2_tof,'pulse time',4*40e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',300e-3,'imaging freq',ImageFreq,'imaging amplitude',1,...
%     'repump shutter delay',2e-3,'imaging_field',dipole_field,'image type','vertical');

% makeImagingSequence_dual_axis(sq,'tof',F2_tof,'tof2',F1_minus_F2_tof,'pulse time',4*40e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',300e-3,'imaging freq',ImageFreq,'imaging amplitude',1,...
%     'repump shutter delay',2e-3,'imaging_field',dipole_field);

sq.find('H-Bridge Helm').before(0.2,0);

setSafeValues(sq);

if nargout == 0
    r = RemoteControl;
    r.upload(sq.compile);
    r.run;
else
    varargout{1} = sq;
end

end