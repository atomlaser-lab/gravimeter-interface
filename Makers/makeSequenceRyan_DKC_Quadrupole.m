function varargout = makeSequenceRyan_DKC_Quadrupole(varargin)   
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
ImageAmp = 0.2;
%% Initialize sequence
sq = initSequence;  %load default values (OLD MOT values are default) 
sq.find('87 imag freq').set(ImageFreq);
sq.find('87 imag amp').set(1);

%% MOT loading
if opt.stage.use_mot
    sq.find("3DMOT").set(0);
    sq.delay(0.5);
    sq.find('2DMOT Freq').set(18);
    sq.find('Push Freq').set(5);
    sq.find('Push amp').set(3.75);
    sq.find('2DMOT').set(1);
    sq.find('3DMOT').set(1);
    sq.find('87 push').set(1);
    % 3D MOT beam settings
    sq.find('3DMOT Freq').set(18);
    sq.find('3DMOT amp').set(1);
    % 3D repump beam settings
    sq.find('87 repump').set(1);
    sq.find('Repump shutter').set(1);
    sq.find('Repump Switch').set(0);
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
    sq.find('87 repump').set(0);
    sq.find('87 repump freq').set(20);
    sq.find('3DMOT freq').set(75);
    sq.delay(Tdepump);
    sq.find('3DMOT').set(0);
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
sq.find('RF atten').set(0);
sq.find('RF Frequency').set(20);
sq.find('Raycus CW').set(0);
sq.find('Raycus TTL').set(0);
sq.find('RedPower CW').set(0);
sq.find('RedPower TTL').set(0);

%% Stern-Gerlach
% sq.delay(20e-3);
% sq.find('CD0 Fast').set(130);
% sq.delay(20e-3);
% sq.find('CD0 Fast').set(0);

%% Take Absorption Image
sq.anchor(timeAtDrop);
sq.camDelay = timeAtDrop - 3;

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