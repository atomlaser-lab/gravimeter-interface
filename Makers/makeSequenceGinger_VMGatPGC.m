function varargout = makeSequenceGinger(varargin)   
%% Parse input arguments
opt = SequenceOptions('load_time',15,'detuning',0,'tof',20e-3,'redpower',2,...
    'keopsys',2);
%keyopsis bec 0.8 0.75
%redpower bec 1.34

% if nargin == 1
%     if ~isa(varargin{1},'SequenceOptions')
%         error('If using only one argument it must of type SequenceOptions');
%     end
%     opt.replace(varargin{1});
% elseif mod(nargin,2) == 0
%     opt.set(varargin{:});
% elseif mod(nargin - 1,2) == 0 && isa(varargin{1},'SequenceOptions')
%     opt.replace(varargin{1});
%     opt.set(varargin{2:end});
% else
%     error('Either supply a single SequenceOptions argument, or supply a set of name/value pairs, or supply a SequenceOptions argument followed by name/value pairs');
% end

% check if SequenceOptions already exists in the workspace
if ~exist('SequenceOptions','class')
    % create a new instance of SequenceOptions if it does not exist
    opt = SequenceOptions('load_time',5,'detuning',0,'tof',17.3e-3,'redpower',2,'keopsys',2);
    %     assignin('caller', 'opt', opt);
else
    % use the existing instance of SequenceOptions
    opt = SequenceOptions();
end

% check input arguments and update options accordingly
if nargin == 0
    if evalin('base', 'exist(''opt'', ''var'')') == 1
        % use the existing instance of SequenceOptions
        opt = evalin('base', 'opt');
    else
        opt = SequenceOptions('load_time', 5, 'detuning', 0, 'tof', 17.3e-3, 'redpower', 2, 'keopsys', 2);
        assignin('base', 'opt', opt);
    end

elseif nargin == 1
    if ~isa(varargin{1},'SequenceOptions')
        error('If using only one argument it must of type SequenceOptions');
    end
    opt.replace(varargin{1});
elseif mod(nargin,2) == 0
    opt.set(varargin{:});
elseif mod(nargin - 1,2) == 0 && isa(varargin{1},'SequenceOptions')
    opt.replace(varargin{1});
    opt.set(varargin{2:end});
    assignin('base', 'opt', opt);

else
    error('Either supply a single SequenceOptions argument, or supply a set of name/value pairs, or supply a SequenceOptions argument followed by name/value pairs');
end

if nargout == 0
    r.make(opt).urun(@Abs_Analysis_Fancy);
    assignin('base', 'opt', opt);
else
    varargout{1} = opt;
end

ImageFreq = opt.detuning*0.6238/6 + 8.3; %assuming imaging field is set to 5V
dipole_field = 5;
imaging_field = dipole_field; %5
ImageAmp = 9; %5 until and including mag load, 9 for high OD samples on horizontal imaging

%make SLM pattern
l=1; %LG charge
phase=0;
% phase=opt.params;
sign=[-1 -1 -1];
f= 79; %79 for charge 1 %focal length
number_of_pulses = 1;
dir = 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\Image Files\512\test\'; %to save images to be sent to the SLM
LG_grating_generator(l,phase,sign,f,dir,number_of_pulses);

%% Initialize sequence
sq = initSequence;  %load default values (OLD MOT values are default) 
sq.find('87 imag freq').set(8.35);
sq.find('87 imag amp').set(8);
%% MOT loading
%
% We use CD channel 0b10 = 2 for loading the MOT 
%
%

sq.find('C6 - N/C').set(1).after(100e-6,0); %SLM trigger

% sq.find('Earth Bias 1').set(7);
% sq.find('Earth Bias 2').set(02.1);
% sq.find('Earth Bias 3').set(7);

sq.find('2DMOT').set(1);
sq.find('3DMOT').set(1);
sq.find('87 push').set(1);
% 3D MOT beam settings
sq.find('3DMOT Freq').set(FtoV('trap',24)); %26 %23 %15.7412
sq.find('3DMOT amp').set(TrapPtoV('trap',1)); %1 %0.77 %1
% 3D repump beam settings
sq.find('87 repump').set(1);
sq.find('Repump Switch').set(0);
sq.find('87 repump freq').set(FtoV('repump',0)); %0
sq.find('87 repump amp').set(TrapPtoV('repump',0.93)); %1
% 3D coil settings
sq.find('H-Bridge Quad').set(1);
sq.find('CD bit 0').set(0); 
sq.find('CD bit 1').set(1);
sq.find('CD2').set(dBtoV('normal',14)); %Coarse control of 3D coils %11
sq.find('CD Fine/Fast').set(dBtoV('fine',8)); %8 fine control of 3D coils
% sq.find('85 repump amp').set(5.5); % 3DMOT bias
% sq.find('MOT bias').set(1); %switching imaging coils on
% sq.find('MOT bias coil').set(1.5); %imaging coils
%Delay for the load time
sq.delay(opt.load_time);

%
% Turn off the 2D MOT and coils as well as the push beam
%
sq.find('2D MOT Coils').before(10e-3,0);
sq.find('2DMOT').before(10e-3,0);
sq.find('87 push').before(10e-3,0);
sq.find('85 push').before(10e-3,0);

sq.find('C6 - N/C').set(1).after(100e-6,0); %setup SLM

%% CMOT sequence
if(1)
%
% Apply a compressed MOT sequence to temporarily increase the density by
% reducing spontaneous emission.  We switch to CD channel 0b00 = 0 because
% it is the fast channel
%

Tcmot = 15e-3; %15e-3
t = 0:1e-3:Tcmot;
%3D Coils
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('CD0 Fast').set(dBtoV('normal',0));
sq.find('CD Fine/Fast').set(dBtoV('fine',12)); 

%Trapping light
sq.find('3DMOT freq').after(t,sq.linramp(t,sq.find('3DMOT freq').values(end),FtoV('trap',46))); %45
sq.find('3DMOT amp').set(TrapPtoV('trap',1)); %1

%Repump
sq.find('87 repump freq').set(FtoV('repump',-8.5)); %-7.5 %-4
sq.find('87 repump amp').set(TrapPtoV('repump',0.05)); %0.08 %0.8

sq.delay(Tcmot);
end
%% PGC sequence
if(1)
%
% Apply polarization gradient cooling to reduce the temperature of the
% atoms.  We use CD channel 0b00 = 0 as it is the fast-switching channel
%

Tpgc = 5e-3; %22e-3
t = 0:1e-3:Tpgc;
sq.find('CD fine/fast').set(dBtoV('fine',13)); %7
sq.find('CD0 Fast').set(dBtoV('normal',0)); %0
sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',72))); %75 %56 %70
sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),TrapPtoV('trap',1))); %1

sq.find('87 repump freq').set(FtoV('repump',-4.8)); %-4.6 %-9 %-2.5 
sq.find('87 repump amp').set(TrapPtoV('repump',0.004)); %0.005 %0.01 %0.05
sq.delay(Tpgc);
end
%% Optical pump atoms into the F = 1 manifold
if(1)
%
% Turn off repump field so that atoms are optically pumped into the F = 1
% manifold.
%
Tdepump = 1e-3;
sq.find('repump switch').set(1); %fiber switch off (it's inverted)
sq.find('87 repump').set(0);
sq.find('87 repump amp').set(0);
sq.find('85 repump').set(0);
sq.find('variable wave plate').set(-2.42); %remove if not doing VMG after PGC
sq.delay(Tdepump);
sq.find('3DMOT').set(0);
end
%% Load into magnetic trap
if(0)
% Load into the magnetic trap at a high gradient.  We switch quickly to a
% low value and then ramp up to the target value
%

Tmagload = 150e-3; %150e-3
t = linspace(0,Tmagload,50);
dBLoad = 110; %110
sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
sq.find('CD Fine/Fast').set(dBtoV('fine',0));
sq.delay(Tmagload); %Tmagload

Toptload = 400e-3; %400e-3
t = linspace(0,Toptload,50);
sq.find('Keopsys MO').set(3.9);
sq.find('Keopsys FA').after(t,sq.minjerk(t,0,DipolePtoV('Keopsys',12))); %5
sq.find('Redpower TTL').set(1);
sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',15))); %15
sq.find('MOT bias').set(1);
sq.find('MOT bias coil').after(t,sq.linramp(t,0,dipole_field));
sq.delay(max(Toptload,Tmagload));

end
%% RF evaporation
if (0)
%
% Remove hot atoms from the sample using RF transitions between the trapped
% |F = 1, m_F = -1> state and the untrapped |F = 1, m_F = 0> state.  All
% frequencies are in MHz
%
rf_start = 20;  
rf_end = 2;
% rf_end = 2; % 20; %1.6 if dipoles are correctly positioned
rf_rate = 2.5; %MHz/s
Tevap = (rf_start - rf_end)/rf_rate;
rf_ramp_type = 'lin';
rf_exp_time_constant = 2;
t = linspace(0,Tevap,100);

sq.find('RF atten').set(1);
if strcmpi(rf_ramp_type,'exp')
    sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
elseif strcmpi(rf_ramp_type,'lin')
    sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
end
sq.delay(Tevap);

sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));

% % % Turn off magnetic trap (use for dipole alignment)
% sq.find('CD0 Fast').set(0);
% sq.find('CD Fine/Fast').set(0);
% sq.delay(30e-3 - opt.tof);
end
%% Ramp down coils
if(0)
%for main coils
Trampcoils = 0.9; %0.5
dB_weak = 0;
% rf_final = 1;
t = linspace(0,Trampcoils,51);
sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));
% sq.find('RF frequency').after(t,sq.linramp(t,sq.find('Rf frequency').values(end),RFtoV(rf_final)));
sq.delay(Trampcoils);

% % for bias coils
% sq.find('85 repump amp').set(0); % 3DMOT bias
% sq.find('MOT bias').set(0); %switching imaging coils on
% sq.find('MOT bias coil').set(0); %imaging coils

% to hold in dipole trap
sq.delay(1-Trampcoils);
end
%% blow away F=2
if(0)
%blow away atoms in F=2
sq.find('87 imag').set(1);
sq.find('variable wave plate').set(-2.42);
sq.delay(1e-3);
sq.find('87 imag').set(0);
end
%% Switch to Helmholtz configuration for state preparation
if (0)
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);

sq.find('H-Bridge Quad').set(0);
sq.delay(50e-6);
sq.find('H-Bridge Helm').set(1);
sq.delay(100e-6);
Tramp = 1;
t = linspace(0,Tramp,51);
sq.find('CD0 Fast').after(t,sq.minjerk(t,0,dBtoV('normal',20)));
%sq.find('MOT bias coil').after(t,sq.minjerk(t,sq.find('MOT bias coil').values(end),0));
sq.delay(Tramp);
end
%% Optical evaporation
if (0)
Tevap = 4; %4 %3
t = linspace(0,Tevap,150);
final_dipole.RP = opt.redpower;
final_dipole.FA = opt.keopsys; 
sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),0.64)); %0.48 %0.42
sq.find('Keopsys FA').after(t,sq.expramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',final_dipole.FA),0.4)); %0.39 %0.42
sq.delay(Tevap);
% sq.delay(1);
end
%% Trigger the DDS
if (1)
sq.ddsTrigDelay = sq.time;
sq.find('DDS TTL').before(10e-3,1).after(10e-3,0);%.after(1e-3,1);
end
%% ARP with DDS
if (0)
sq.dds(1).set(38,0,0);
sq.dds(2).set(110,0,0);
sq.delay(10e-6);
sq.find('RF Switch').set(1);
Tarp = 80e-3; %80e-3
t = linspace(0,Tarp,501);
%df = 38.35 31/10/22
df = 38.3 + 1*sq.linramp(t,-0.5,0.5);
w = Tarp/10;
% amp = 0.075*sech((t - Tarp/2)/w).^2;
% amp = 1.3e-3*sech((t - Tarp/2)/w).^2;
amp = 1.5e-3*sech((t - Tarp/2)/w).^2;
% amp = 0.0043*ones(size(t));
sq.dds(1).after(t,df,amp,0);
sq.dds(2).after(t,110,0,0);
sq.delay(Tarp);
sq.dds(1).set(110,0,0);
sq.dds(2).set(110,0,0);
sq.find('RF Switch').set(0);
end
%% RF Pi Pulse from |1,-1> to |1,0>
if(1)
pulse_freq = 38.336; %38.263 was quoted to work on 16/9/22
sq.dds(1).set(pulse_freq,0,0);
sq.dds(2).set(110,0,0);
sq.delay(20e-6);
sq.find('RF Switch').set(1);
sq.dds(1).set(pulse_freq,0.075*1,0);
sq.dds(2).set(110,0,0);
sq.delay(10e-6);
sq.dds(1).set(pulse_freq,0,0);
sq.dds(2).set(110,0,0);
sq.find('RF Switch').set(0);
end
%% Drop atoms
% sq.delay(1);
% sq.find('Earth Bias 1').set(0);
% sq.find('Earth Bias 2').set(0);
% sq.find('Earth Bias 3').set(0);

timeAtDrop = sq.time;
% sq.find('Probe').set(1).after(1e-3,0);
%
% This trigger delay is necessary because the DDS instructions start when
% the DDS trigger occurs
%
% sq.ddsTrigDelay = timeAtDrop;
% sq.find('DDS TTL').before(10e-3,1).after(10e-3,0).after(1e-3,1);
%
% Set all other channels to 0

sq.find('2D MOT Coils').set(0);
sq.find('3DMOT').set(0);
sq.find('87 repump').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(dBtoV('normal',0)); %REMEMBER TO TURN BACK TO 2 WHEN USING RAMAN??
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

%marker for bottom of cell at 45ms from drop (22ms for horizontal imaging)
% sq.find('C6 - N/C').set(0).after(45e-3,1).after(10e-6,0);

% sq.delay(100e-6);
% sq.find('H-bridge helm').set(0);
% sq.delay(50e-6);
% sq.find('H-bridge quad').set(1);

%% VMG!
if(1)
sq.anchor(timeAtDrop);

raman_delay_drop = 5e-3; 
% raman_delay_drop = opt.params;
interrogation_time = 10e-6;
% interrogation_time = opt.params;
Traman_pi = 100e-6;
% Traman_pi = opt.params;
% F2_tof = 2e-3; %12ms for horizontal, 6ms for vertical
F2_tof = raman_delay_drop + Traman_pi + 3e-4;
% F2_tof = raman_delay_drop + Traman_pi/2 + interrogation_time + Traman_pi/2 + 10e-6;
F1_tof = F2_tof + 3e-3;
F1_minus_F2_tof = F1_tof - F2_tof; %needs to be 12ms minimum for full frame
% AI_phi = opt.params;
AI_phi = 0;
sq.anchor(timeAtDrop + raman_delay_drop);
Power_raman_G = 1; %1 %sideband
Power_raman_LG = 0.1; %0.1 %carrier 
% Power_raman_LG = opt.params;
% Delta_raman = opt.params;
Delta_raman = 19.832;
 
% state prep
sq.dds(1).set(110+Delta_raman/4,Power_raman_G,0); 
sq.dds(2).set(110-Delta_raman/4,Power_raman_LG,0); 
sq.delay(Traman_pi);
sq.dds(1).set(110,0,0);
sq.dds(2).set(110,0,0);
sq.delay(10e-6);
sq.find('C6 - N/C').set(1).after(100e-6,0); %SLM trigger

% %blow away F=1 atoms and give SLM enough time to be triggered again
% sq.find('85 imag').set(1);
% sq.find('85 imag freq').set(0); %0
% sq.find('85 imag amp').set(8); %1
% sq.delay(1e-3);
% sq.find('85 imag').set(0); %1

% % BS1
% sq.dds(1).set(110+Delta_raman/4,Power_raman_G,0); %sideband
% sq.dds(2).set(110-Delta_raman/4,Power_raman_LG,0); %carrier
% sq.delay(Traman_pi/2);
% % sq.delay(Traman_pi);
% sq.dds(1).set(110,0,0);
% sq.dds(2).set(110,0,0);
% sq.delay(10e-6);
% sq.find('C6 - N/C').set(1).after(100e-6,0); %SLM trigger
% 
% sq.delay(interrogation_time);
% 
% % BS2
% sq.dds(1).set(110+Delta_raman/4,Power_raman_G,0); %sideband
% sq.dds(2).set(110-Delta_raman/4,Power_raman_LG,AI_phi); %carrier
% sq.delay(Traman_pi/2);
% % sq.delay(Traman_pi);
% sq.dds(1).set(110,0,0);
% sq.dds(2).set(110,0,0);

% blow away atoms in F=2
% sq.delay(0.5e-3);
% sq.find('87 imag').set(1);
% sq.delay(1e-3);
% sq.find('87 imag').set(0);
end
%% S -G Field (SG pulse)
if(0)

% HEY RYAN! MAKE SURE THIS IS UNCOMMENTED TO SWITCH BACK TO ANTI-HELMHOLTZ CONFIGURATION FOR STERN-GERLACH
sq.delay(10e-6);
sq.find('CD0 Fast').set(0);
sq.delay(100e-6);
sq.find('H-bridge helm').set(0);
sq.delay(50e-6);
sq.find('H-bridge quad').set(1);

sq.anchor(timeAtDrop);
sq.delay(1e-3); %raman_delay_drop+1e-3
SG_Pulse_time = 8e-3;
SG_Amp = 25; %25
t = linspace(0,SG_Pulse_time,50);
sq.find('CD0 Fast').after(t,sq.minjerk(t,0,dBtoV('normal',SG_Amp)));
sq.delay(SG_Pulse_time);
sq.find('CD0 Fast').after(t,sq.minjerk(t,sq.find('CD0 Fast').values(end),0));

% F2_tof = 3e-3;
% F1_tof = 15e-3;
end
%% Take Absorption Image

sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 2;

% makeImagingSequence(sq,'tof',opt.tof,'pulse time',40e-6,'repump delay',100e-6,...
%     'repump time',200e-6,'cam time',5e-6,'cycle time',100e-3,...
%     'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
%     'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal'); 
 
makeImagingSequence_4_Images(sq,'tof',F2_tof,'tof2',F1_minus_F2_tof,'pulse time',40e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',50e-6,'cycle time',100e-3,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
    'fibre switch delay',1e-3,'imaging_field',imaging_field,'image type','horizontal');

% turn off the dipoles
sq.find('Redpower CW').set(0);
sq.find('Redpower TTL').after(100e-6,0);
sq.find('Keopsys FA').set(0);
sq.find('Keopsys MO').after(100e-6,0);
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));
sq.find('3DMOT').set(0);
sq.find('87 repump amp').set(0);
sq.find('CD0 Fast').set(0);
% sq.find('C6 - N/C').set(1).after(100e-6,0).after(10e-3,1).after(100e-6,0);
% sq.delay(100e-6);
% sq.find('H-bridge helm').set(0);
% sq.delay(50e-6);
% sq.find('H-bridge quad').set(1);

sq.waitFromLatest(0.25);
setSafeValues(sq);
% sq.delay(5);
 
if nargout == 0
    r = RemoteControl;
    r.upload(sq.compile);
    r.run;
else
    varargout{1} = sq;
end

end