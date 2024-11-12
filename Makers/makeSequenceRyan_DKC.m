function varargout = makeSequenceRyan_DKC(varargin)   
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

% ImageFreq = opt.detuning*0.6238/6 + 8.3;  %For a 0.5 V bias field during imaging
ImageFreq = FtoV('image',opt.detuning + 0*3.1); %at 3 V imaging field, 3.1 MHz detuning, 2.73 MHz at 10 V
% ImageFreq = opt.params(1);
dipole_field = 5;
% dipole_field = opt.params(1);
imaging_field = 0.5;
ImageAmp = 7;
%% Initialize sequence
sq = initSequence;  %load default values (OLD MOT values are default) 
sq.find('87 imag freq').set(ImageFreq);
sq.find('87 imag amp').set(8);
% sq.find('Variable wave plate').set(-5); %Set for ND imaging
%% MOT loading
%
% We use CD channel 0b10 = 2 for loading the MOT 
%
sq.find('2DMOT').set(1);
sq.find('3DMOT').set(1);
sq.find('87 push').set(1);
% 3D MOT beam settings
sq.find('3DMOT Freq').set(FtoV('trap',24));
sq.find('3DMOT amp').set(TrapPtoV('trap',1));
% 3D repump beam settings
sq.find('87 repump').set(1);
sq.find('Repump Switch').set(0);
sq.find('87 repump freq').set(FtoV('repump',0));
sq.find('87 repump amp').set(TrapPtoV('repump',1));
% 3D coil settings
sq.find('H-Bridge Quad').set(1);
sq.find('CD bit 0').set(0); 
sq.find('CD bit 1').set(1);
sq.find('CD2').set(dBtoV('normal',24)); %Coarse control of 3D coils
sq.find('CD Fine/Fast').set(dBtoV('fine',0)); % fine control of 3D coils
% Bias coil settings
sq.find('Earth Bias 1').set(3); %6
sq.find('Earth Bias 2').set(5); %4
sq.find('Earth Bias 3').set(0); %2
%Delay for the load time
sq.delay(opt.load_time);
%
% Turn off the 2D MOT and coils as well as the push beam
%
sq.find('2D MOT Coils').before(10e-3,0);
sq.find('2DMOT').before(10e-3,0);
sq.find('87 push').before(10e-3,0);
sq.find('85 push').before(10e-3,0);

%% CMOT sequence
%
% Apply a compressed MOT sequence to temporarily increase the density by
% reducing spontaneous emission.  We switch to CD channel 0b00 = 0 because
% it is the fast channel
%
Tcmot = 30e-3;
t = 0:1e-3:Tcmot;
%3D Coils
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('CD0 Fast').set(dBtoV('normal',0));
sq.find('CD Fine/Fast').set(dBtoV('fine',10)); 
%Trapping light
sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',53.4)));
sq.find('3DMOT amp').set(TrapPtoV('trap',1));
%Repump
sq.find('87 repump freq').set(FtoV('repump',12)); %-7
sq.find('87 repump amp').set(TrapPtoV('repump',0.05)); %0.05

sq.delay(Tcmot);

%% PGC sequence
%
% Apply polarization gradient cooling to reduce the temperature of the
% atoms.  We use CD channel 0b00 = 0 as it is the fast-switching channel
%
Tpgc = 6e-3;
t = 0:0.5e-3:Tpgc;
% t = linspace(0,Tpgc,26);
sq.find('Earth Bias 1').set(1.2); %Probably N/S
sq.find('Earth Bias 2').set(6.5);%Pr0bably U/D
sq.find('Earth Bias 3').set(0.05); % Probably E/W
sq.find('CD fine/fast').set(dBtoV('fine',0));
sq.find('CD0 Fast').set(dBtoV('normal',0));
sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',79.8))); 
sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),TrapPtoV('trap',0.8))); %0.5

sq.find('87 repump freq').set(FtoV('repump',11));%-4.8
sq.find('87 repump amp').set(TrapPtoV('repump',0.05));%0.004

sq.delay(Tpgc);

%% Optical pump atoms into the F = 1 manifold
%
% Turn off repump field so that atoms are optically pumped into the F = 1
% manifold.
%
Tdepump = 1e-3;
sq.find('repump switch').set(1); %fiber switch off (it's inverted)
sq.find('87 repump').set(0);
sq.find('87 repump amp').set(0);
% sq.find('85 repump').set(0);
% sq.find('85 repump amp').set(0);
sq.delay(Tdepump);
sq.find('3DMOT').set(0);

%% Load into magnetic trap
%
% Load into the magnetic trap at a high gradient.  We switch quickly to a
% low value and then ramp up to the target value
%
Tmagload = 150e-3;
t = 0:10e-3:Tmagload;
sq.find('Earth Bias 1').set(0); %Probably N/S
sq.find('Earth Bias 2').set(0);%Pr0bably U/D
sq.find('Earth Bias 3').set(0); % Probably E/W
dBLoad = 110;
sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
sq.find('CD Fine/Fast').set(dBtoV('fine',0));
%  sq.delay(Tmagload);

Toptload = 400e-3;
t = 0:20e-3:Toptload;
sq.find('Raycus TTL').set(1);
sq.find('Raycus CW').after(t,sq.minjerk(t,0,DipolePtoV('raycus',10)));
sq.find('Redpower TTL').set(1);
sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',15)));
sq.find('MOT bias').set(1);
sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));
sq.delay(max(Toptload,Tmagload));

%% RF evaporation
%
% Remove hot atoms from the sample using RF transitions between the trapped
% |F = 1, m_F = -1> state and the untrapped |F = 1, m_F = 0> state.  All
% frequencies are in MHz
%
rf_start = 16;
rf_end = 2;
rf_rate = 3.5;    %MHz/s
Tevap = (rf_start - rf_end)/rf_rate;
rf_ramp_type = 'lin';
rf_exp_time_constant = 2;
t = linspace(0,Tevap,50);

sq.find('RF atten').set(1);
if strcmpi(rf_ramp_type,'exp')
    sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
elseif strcmpi(rf_ramp_type,'lin')
    sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
end
sq.delay(Tevap);

sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));

%% Turn off magnetic trap
% Trampcoils = 0.8;
% dB_weak = 0;
% t = linspace(0,Trampcoils,51);
% sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));
% sq.delay(Trampcoils);

sq.find('CD0 Fast').set(0);
sq.find('CD Fine/Fast').set(0);
% sq.delay(15e-3 - opt.tof);

%% Optical evaporation
Tevap = 2.5 + 0.0;
TC = 0.5;
t = linspace(0,Tevap,150);
final_dipole.RP = opt.redpower;
final_dipole.FA = opt.raycus;
sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),TC));
sq.find('Raycus CW').after(t,sq.expramp(t,sq.find('Raycus CW').values(end),DipolePtoV('raycus',final_dipole.FA),TC));
sq.delay(Tevap);

if ~isempty(opt.dkc.dipoles)
    sq.find('Raycus CW').set(DipolePtoV('raycus',opt.dkc.dipoles(1)));
    sq.find('RedPower CW').set(DipolePtoV('redpower',opt.dkc.dipoles(2)));
    sq.delay(500e-3);
end

%% Drop atoms
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
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));
% sq.find('Variable wave plate').set(-3.7);   %This value switches to absorption imaging

%% Apply guiding
if opt.dkc.is_guide()
    sq.find('dkc power').set(DipolePtoV('dkc',opt.dkc.power)).before(10e-3,0.3);
    sq.delay(opt.dkc.duration);
    sq.find('dkc power').set(0);
end

%% Take Absorption Image
sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 2;

makeImagingSequence(sq,'tof',opt.tof,'pulse time',80e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',5e-6,'cycle time',50e-3,...
    'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
    'fibre switch delay',1e-3,'imaging_field',dipole_field,'image type','horizontal');

% turn off the dipoles
sq.find('Redpower CW').set(0);
sq.find('Redpower TTL').after(100e-6,0);
sq.find('Raycus CW').set(0);
sq.find('Raycus TTL').after(100e-6,0);
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));
sq.find('3DMOT').set(0);
sq.find('87 repump amp').set(0);

%% Apply DKC
if opt.dkc.is_kick()
    sq.anchor(timeAtDrop);
    sq.find('dkc power').before(10e-3,0.3);
    sq.delay(opt.dkc.delay);
    sq.find('dkc power').set(DipolePtoV('dkc',opt.dkc.power));
    sq.delay(opt.dkc.duration);
    sq.find('dkc power').set(0);
    sq.find('dkc ttl').set(0);
end

%% Reproduce guiding/DKC effect for image correction
% sq.waitFromLatest(0.1);
% sq.find('dkc power').set(0.35);
% sq.delay(0.5);
% sq.find('dkc power').set(DipolePtoV('dkc',opt.params(1)));
% sq.delay(20e-3);
% sq.find('dkc power').set(0);
% sq.delay(opt.tof - 20e-3);
% sq.find('vertical cam trig').set(1).after(10e-6,0);
% sq.delay(15e-6);
% sq.find('87 imag').set(1);
% sq.delay(40e-6);
% sq.find('87 imag').set(0);


%% Finish
if opt.dkc.is_test()
    sq.find('dkc power').set(0.36);
    sq.waitFromLatest(0.25);
    sq.find('vertical cam trig').set(1);
    sq.delay(10e-6);
    sq.find('dkc power').set(0);
    sq.find('vertical cam trig').set(0);
    sq.delay(100e-6);
    sq.find('dkc ttl').set(0);
end

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