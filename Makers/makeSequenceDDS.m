function varargout = makeSequenceDDS(varargin)
opt = parse_maker_variable_argument_list(varargin{:});
sq = initSequence;  %load default values (OLD MOT values are default)
sq.find('2DMOT').set(0);
sq.find('3DMOT').set(0);
sq.find('87 push').set(0);
sq.find('87 repump').set(0);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('Vertical MOT Mirror').set(1);
sq.delay(1);

sq.ddsTrigDelay = sq.time;
sq.find('DDS Trigger').before(10e-3,1).after(10e-3,0);%.after(1e-3,1);

sq.delay(100e-3);

timeAtDrop = sq.time;

sq.anchor(timeAtDrop);

% state prep
% sq.dds(1).
Tarp = 10;
t = linspace(0,Tarp,100);
f = linspace(-6,6,100);
amp = sq.linramp(t,0,1);
sq.dds(1).after(t,110,0,0);
sq.dds(2).after(t,110 + f,1,0);
% sq.delay(Traman_pi);
sq.dds(1).set(110,0,0);
sq.dds(2).set(110,0,0);

sq.delay(10e-6);
makeImagingSequence_4_Images(sq,'tof',10e-3,'tof2',3e-3,'pulse time',4*40e-6,'repump delay',100e-6,...
    'repump time',200e-6,'cam time',5e-6,'cycle time',300e-3,'imaging freq',0,'imaging amplitude',1,...
    'repump shutter delay',2e-3,'imaging_field',1,'image type','vertical');

setSafeValues(sq);

if nargout == 0
    r = RemoteControl;
    r.upload(sq.compile);
    r.run;
else
    varargout{1} = sq;
end

end