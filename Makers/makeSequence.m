function varargout = makeSequence(varargin)   
%% Parse input arguments
opt = parse_maker_variable_argument_list(varargin{:});

% ImageFreq = opt.detuning;
% ImageFreq = opt.detuning + 2.3; %For laser cooling stages
ImageFreq = opt.detuning + 0.5; %Low intensity after dipole evaporation
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
    sq.find('H-Bridge Quad').set(1);
    sq.find('H-Bridge Helm').set(0);
    sq.find('CD bit 0').set(0);
    sq.find('CD bit 1').set(0);
    sq.find('CD0 Fast').set(14); %Coarse control of 3D coils
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
    dBmax = 110;
    dBLoad = 55;
    sq.find('CD0 Fast').after(t,sq.linramp(t,dBLoad,dBmax));
    sq.find('CD Fine/Fast').set(0);
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
    rf_start = 16; %16
    rf_end = 0.75;
%     rf_end = 4;
    rf_rate = 3;    %MHz/s 3
    Tevap = (rf_start - rf_end)/rf_rate;
    t = linspace(0,Tevap,50);
    
    sq.find('RF switch').set(1); %NOTE: FG/DDS is set to 1 in initSequence, so FG is the default RF source
    sq.find('RF frequency').set(rf_start);
    sq.delay(0.25);
    sq.find('RF frequency').after(t,sq.linramp(t,rf_start,rf_end));
    sq.delay(Tevap);
    
    sq.find('RF switch').set(0);
    sq.find('RF Frequency').set(20);
end

%% Test loading atoms into magnetic trap
% sq.find('CD0 Fast').set(0);
% sq.delay(35e-3 - opt.tof);
% sq.find('Raycus CW').set(0);
% sq.find('Raycus TTL').set(0);
% sq.find('RedPower CW').set(0);
% sq.find('RedPower TTL').set(0);


%% Load into dipole trap
if opt.stage.use_dipoles
    Trampcoils = 0.3;
    dB_weak = 0;
    t = linspace(0,Trampcoils,51);
    sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dB_weak));
    sq.find('MOT bias coil').after(t,sq.linramp(t,sq.find('MOT bias coil').values(end),dipole_field));
    sq.delay(Trampcoils);
end

%% Optical evaporation
if opt.stage.use_evap_dipoles
    Tevap = 5;
    TC = 0.75;
    t = linspace(0,Tevap,51);
    sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),opt.redpower,TC));
    sq.find('Raycus CW').after(t,sq.expramp(t,sq.find('Raycus CW').values(end),opt.raycus,TC));
    sq.delay(Tevap);
end

%%
% sq.delay(0.1);

%% Drop atoms
timeAtDrop = sq.time;
sq.find('2D MOT Coils').set(0);
sq.find('3DMOT').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(0);
sq.find('CD2').set(0);
sq.find('CD Fine/Fast').set(0);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
% sq.find('H-Bridge Quad').set(0);
% sq.find('H-Bridge Helm').set(0);
sq.find('RF switch').set(0);
sq.find('RF Frequency').set(20);
sq.find('Raycus CW').set(0);
sq.find('Raycus TTL').set(0);
sq.find('RedPower CW').set(0);
sq.find('RedPower TTL').set(0);

%% Stern-Gerlach
% sq.delay(10e-3);
% sq.find('CD0 Fast').set(100);
% sq.delay(10e-3);
% sq.find('CD0 Fast').set(0);

%% Take Absorption Image
sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 3;

makeImagingSequence(sq,'tof',opt.tof,'pulse time',4*40e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
    'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
    'repump shutter delay',2e-3,'imaging_field',dipole_field,'image type','horizontal');


setSafeValues(sq);
 
if nargout == 0
    r = RemoteControl;
    r.upload(sq.compile);
    r.run;
else
    varargout{1} = sq;
end

end