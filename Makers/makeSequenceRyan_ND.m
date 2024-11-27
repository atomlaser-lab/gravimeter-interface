function varargout = makeSequenceRyan_ND(varargin)   
%% Parse input arguments
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

ImageFreq = opt.detuning + 1.0784;
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

if opt.nd.enable_fb_laser
    sq.find('Feedback Laser TTL').set(1);
    sq.find('Feedback laser power').set(200e-3);
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
    sq.find('2D MOT Coils').before(10e-3,0);
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
%     sq.find('87 repump amp').set(0).after(5e-3,1);
    sq.find('87 repump freq').set(0);
    sq.find('3DMOT freq').set(75);
%     sq.find('MOT bias coil').before(1e-3,3);
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
%     sq.find('Bias E/W').set(0);
%     sq.find('Bias N/S').set(0);
%     sq.find('Bias U/D').set(0);
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
        sq.find('Raycus CW').after(t,sq.linramp(t,sq.find('Raycus CW').values(end),4));
        sq.find('Redpower CW').after(t,sq.linramp(t,sq.find('Redpower CW').values(end),12));
        if opt.nd.enable_fb_laser
            sq.find('Feedback Laser TTL').set(1);
            sq.find('Feedback laser power').after(t,sq.linramp(t,sq.find('Feedback laser power').values(end),opt.nd.fb_laser_power));
        end
        sq.delay(max(Toptload - Tmagload,0));
    end

    if ~opt.stage.use_evap_mag
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
    rf_start = 16;
    rf_end = 0.75;
    rf_rate = 3;    %MHz/s
    Tevap = (rf_start - rf_end)/rf_rate;
    t = linspace(0,Tevap,50);
    
    sq.find('RF atten').set(1);
    sq.find('RF frequency').set(rf_start);
    sq.delay(0.25);
    sq.find('RF frequency').after(t,sq.linramp(t,rf_start,rf_end));
    sq.delay(Tevap);
    
    sq.find('RF atten').set(0);
    sq.find('RF Frequency').set(20);
end

%% Take dummy images for NDI
%
% This takes 2 images for NDI, where the first is a "dummy" image to get
% the camera to properly time its acquisition, and the second is a
% reference for NDI.  These images are taken 3 seconds before evaporation
% ends.  Sequence time is re-anchored to the time at which evaporation ends
%
if opt.stage.use_dipoles && opt.nd.enable_ndi && opt.nd.ref_images > 0
    time_at_evap_end = sq.time;
    sq.anchor(sq.time - 3);
    sq.camDelay = sq.time - 2;
    makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_time,'cycle time',500e-3,...
        'imaging amplitude',opt.nd.pulse_power,'num_images',opt.nd.ref_images,...
        'pulse delay',opt.nd.pulse_delay);
    sq.anchor(time_at_evap_end);
end

%% Test loading atoms into magnetic trap
% sq.find('CD0 Fast').set(0);
% sq.delay(30e-3 - opt.tof);
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

%% Non-destructive imaging/feedback
if opt.nd.enable_ndi
%     t = 0:5e-3:250e-3;
%     sq.find('Raycus CW').after(t,sq.find('Raycus CW').values(end) + sq.linramp(t,0,0.2));
%     sq.find('RedPower CW').after(t,sq.find('Redpower CW').values(end) + sq.linramp(t,0,0.2));
    sq.delay(250e-3);
    
%     quad_amp = 0.25;%*opt.params(1);
%     quadrupole_freq = opt.params(1);
%     num_cycles = 10;
%     dt = 0.1/quadrupole_freq;
%     T = num_cycles/quadrupole_freq;
%     t = 0:dt:T;
%     ch = sq.find('Feedback laser power');
% %     ch = sq.find('RedPower CW');
%     quad_driving_signal = ch.values(end) + quad_amp*sin(2*pi*quadrupole_freq*t);
%     ch.after(t,quad_driving_signal);
%     sq.delay(T);

    sq.find('Feedback laser power').after(50e-3,1);

    makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_time,'cycle time',opt.nd.cycle_time,...
        'imaging amplitude',opt.nd.pulse_power,'num_images',opt.nd.num_images(1),'pulse delay',opt.nd.pulse_delay);

%     sq.delay(1);
%     makeNDImagingSequence(sq,'pulse time',opt.nd.pulse_time,'cam time',opt.nd.pulse_time,'cycle time',opt.nd.cycle_time,...
%         'imaging amplitude',opt.nd.pulse_power,'num_images',opt.nd.num_images(2),'pulse delay',opt.nd.pulse_delay);
end

%%
% t = 0:5e-3:250e-3;
% sq.find('Raycus CW').after(t,sq.find('Raycus CW').values(end) + sq.linramp(t,0,0.2));
% sq.find('RedPower CW').after(t,sq.find('Redpower CW').values(end) + sq.linramp(t,0,0.2));
% sq.delay(250e-3);
% 
% quad_amp = 0.25*opt.params(1);
% quadrupole_freq = 33;
% % quadrupole_freq = opt.params(1);
% num_cycles = 10;
% dt = 0.1/quadrupole_freq;
% T = num_cycles/quadrupole_freq;
% t = 0:dt:T;
% ch = sq.find('Feedback laser power');
% % ch = sq.find('RedPower CW');
% quad_driving_signal = ch.values(end) + quad_amp*sin(2*pi*quadrupole_freq*t);
% ch.after(t,quad_driving_signal);
% sq.delay(T);
% 
% sq.delay(500e-3);

% sq.delay(0.5);
% sq.find('CD0 Fast').set(0);
% sq.find('Feedback Laser TTL').before(1.5,1);
% sq.find('Feedback laser power').before(1.5,10);
% sq.delay(30e-3 - opt.tof);
% sq.find('Feedback Laser TTL').set(0);
% % sq.delay(1);

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
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(20);
sq.find('Raycus CW').set(0);
sq.find('Raycus TTL').set(0);
sq.find('RedPower CW').set(0);
sq.find('RedPower TTL').set(0);
sq.find('Feedback laser power').set(0);
sq.find('Feedback Laser TTL').set(0);

%% Stern-Gerlach
sq.delay(15e-3);
sq.find('CD0 Fast').set(75);
sq.delay(5e-3);
sq.find('CD0 Fast').set(0);

%% Take Absorption Image
sq.anchor(timeAtDrop);
if opt.nd.ref_images == 0 || opt.nd.enable_ndi == 0
    sq.camDelay = timeAtDrop - 3;
end
sq.find('ND imag amp').set(1);
makeImagingSequence(sq,'tof',opt.tof,'pulse time',40e-6,'repump delay',100e-6,...
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