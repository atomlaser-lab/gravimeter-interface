function varargout = makeSequenceRyan(varargin)   
%% Parse input arguments
opt = SequenceOptions('load_time',7.5,'detuning',0,'tof',20e-3,'redpower',2,...
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

ImageFreq = FtoV('image',opt.detuning);
ImageAmp = 5.0;
%% Initialize sequence
sq = initSequence;  %load default values (OLD MOT values are default) 
%% MOT loading
%
% We use CD channel 0b10 = 2 for loading the MOT 
%
sq.find('2DMOT').set(1);
sq.find('3DMOT').set(1);
sq.find('87 push').set(1);
% 3D MOT beam settings
sq.find('3DMOT Freq').set(FtoV('trap',20));
sq.find('3DMOT amp').set(TrapPtoV('trap',1));
% 3D repump beam settings
sq.find('87 repump').set(1);
sq.find('Repump Switch').set(0);
sq.find('87 repump freq').set(FtoV('repump',55));
sq.find('87 repump amp').set(TrapPtoV('repump',1));
% 3D coil settings
sq.find('H-Bridge Quad').set(1);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(1);
sq.find('CD2').set(dBtoV('normal',11));
sq.find('CD Fine/Fast').set(dBtoV('fine',3.44));
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
%3D Coils
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('CD Fine/Fast').set(dBtoV('fine',8.8)); 
%Trapping light
sq.find('3DMOT freq').set(FtoV('trap',30)); 
sq.find('3DMOT amp').set(TrapPtoV('trap',0.8));
%Repump
sq.find('87 repump freq').set(FtoV('repump',55)); 
sq.find('87 repump amp').set(TrapPtoV('repump',0.1)); 

sq.delay(Tcmot);

%% PGC sequence
%
% Apply polarization gradient cooling to reduce the temperature of the
% atoms.  We use CD channel 0b00 = 0 as it is the fast-switching channel
%
Tpgc = 5e-3;
t = linspace(0,Tpgc,100);
% sq.find('CD fine/fast').set(dBtoV('fine',0));
sq.find('CD fine/fast').after(t,sq.minjerk(t,sq.find('CD fine/fast').values(end),dBtoV('fine',2)));
sq.find('3DMOT freq').after(t,sq.minjerk(t,sq.find('3DMOT freq').values(end),FtoV('trap',60)));
sq.find('3DMOT amp').after(t,sq.minjerk(t,sq.find('3DMOT amp').values(end),TrapPtoV('trap',0.5)));

sq.find('87 repump freq').set(FtoV('repump',55)); 
sq.find('87 repump amp').set(TrapPtoV('repump',1)); 

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
sq.find('85 repump').set(0);
sq.find('85 repump amp').set(0);
sq.delay(Tdepump);
sq.find('3DMOT').set(0);

%% Load into magnetic trap
%
% Load into the magnetic trap at a high gradient.  We switch quickly to a
% low value and then ramp up to the target value
%
Tmagload = 150e-3;
t = linspace(0,Tmagload,100);
dBLoad = 110;
sq.find('CD0 Fast').after(t,sq.linramp(t,dBtoV('normal',dBLoad/2),dBtoV('normal',dBLoad)));
% sq.find('CD0 Fast').set(dBtoV('normal',50));
sq.find('CD Fine/Fast').set(dBtoV('fine',0));
% sq.delay(250e-3);

Toptload = 400e-3;
t = linspace(0,Toptload,100);
% sq.find('Keopsys MO').set(3.9);
% sq.find('Keopsys FA').after(t,sq.minjerk(t,0,DipolePtoV('Keopsys',11)));
% sq.find('Redpower TTL').set(1);
% sq.find('Redpower CW').after(t,sq.minjerk(t,0,DipolePtoV('RedPower',17.5)));
sq.delay(Toptload);

%% RF evaporation
%
% Remove hot atoms from the sample using RF transitions between the trapped
% |F = 1, m_F = -1> state and the untrapped |F = 1, m_F = 0> state.  All
% frequencies are in MHz

% rf_start = 15;
% rf_end = 7;
% rf_rate = 3;    %MHz/s
% Tevap = (rf_start - rf_end)/rf_rate;
% rf_ramp_type = 'lin';
% rf_exp_time_constant = 2;
% t = linspace(0,Tevap,200);
% 
% sq.find('RF atten').set(1);
% if strcmpi(rf_ramp_type,'exp')
%     sq.find('RF frequency').after(t,sq.expramp(t,RFtoV(rf_start),RFtoV(rf_end),rf_exp_time_constant)); %ramp rf frequency from 4 to -2.667
% elseif strcmpi(rf_ramp_type,'lin')
%     sq.find('RF frequency').after(t,sq.linramp(t,RFtoV(rf_start),RFtoV(rf_end)));
% end
% sq.delay(Tevap);

%% Weaken magnetic trap, continue evaporation
% Trampcoils = 0.5;
% dB_weak = 8;
% FinalRFKnife = 15;
% t = linspace(0,Trampcoils,100);
% sq.find('CD0 Fast').after(t,sq.linramp(t,sq.find('CD0 Fast').values(end),dBtoV('normal',dB_weak)));
% sq.find('RF Frequency').after(t,sq.linramp(t,sq.find('RF Frequency').values(end),RFtoV(FinalRFKnife)));
% sq.delay(Trampcoils);
% sq.find('RF atten').set(0);
% sq.find('RF Frequency').set(RFtoV(20));

% Turn off magnetic trap
% sq.find('CD0 Fast').set(0);
% sq.find('CD Fine/Fast').set(0);
% sq.delay(30e-3 - opt.tof);

%% Optical evaporation
% Tevap = 4;
% t = linspace(0,Tevap,200);
% final_dipole.RP = opt.redpower;
% final_dipole.FA = opt.keopsys;
% % final_dipole.RP = 1.25;
% % final_dipole.FA = 1.25;
% sq.find('RedPower CW').after(t,sq.expramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',final_dipole.RP),0.4));
% sq.find('Keopsys FA').after(t,sq.expramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',final_dipole.FA),0.4));
% sq.delay(Tevap);
% 
% T2 = 1;
% t = linspace(0,T2,100);
% sq.find('RedPower CW').after(t,sq.linramp(t,sq.find('RedPower CW').values(end),DipolePtoV('redpower',2)));
% sq.find('Keopsys FA').after(t,sq.linramp(t,sq.find('Keopsys FA').values(end),DipolePtoV('keopsys',2)));
% sq.delay(T2);

%% Drop atoms
sq.delay(1);
timeAtDrop = sq.time;
sq.find('CD0 Fast').set(0);
sq.find('CD Fine/Fast').set(0);
sq.find('Redpower CW').set(0);
sq.find('Redpower TTL').after(100e-6,0);
sq.find('Keopsys FA').set(0);
sq.find('Keopsys MO').after(100e-6,0);
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));

%% Take Absorption Image
sq.anchor(timeAtDrop);
makeImagingSequence(sq,'tof',opt.tof,'pulse time',0.1e-3,'repump delay',100e-6,...
    'repump time',100e-6,'cam time',500e-6,'cycle time',40e-3,...
    'manifold',1,'imaging freq',ImageFreq,'imaging amplitude',ImageAmp,...
    'fibre switch delay',1e-3);

% turn off the dipoles
sq.find('Redpower CW').set(0);
sq.find('Redpower TTL').after(100e-6,0);
sq.find('Keopsys FA').set(0);
sq.find('Keopsys MO').after(100e-6,0);
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(RFtoV(20));
 
if nargout == 0
    r = RemoteControl;
    r.upload(sq.compile);
    %to plot
    sq.plot;
    r.run;
else
    varargout{1} = sq;
end

end