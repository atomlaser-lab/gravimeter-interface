function varargout = makeSequenceDDS(varargin)
%% Parse input arguments
opt = SequenceOptions('load_time',15,'detuning',0,'tof',20e-3,'redpower',2,...
    'raycus',2);

% check if SequenceOptions already exists in the workspace
if ~exist('SequenceOptions','class')
    % create a new instance of SequenceOptions if it does not exist
    opt = SequenceOptions('load_time',5,'detuning',0,'tof',17.3e-3,'redpower',2,'raycus',2);
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
        opt = SequenceOptions('load_time', 5, 'detuning', 0, 'tof', 17.3e-3, 'redpower', 2, 'raycus', 2);
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

% Initialize sequence
sq = initSequence;  %load default values (OLD MOT values are default)

sq.find('2DMOT').set(0);
sq.find('3DMOT').set(0);
sq.find('87 push').set(0);
sq.find('87 repump').set(0);
sq.find('Repump Switch').set(1);
sq.find('H-Bridge Quad').set(0);
sq.find('CD bit 0').set(0);
sq.find('CD bit 1').set(0);
sq.find('Earth Bias 1').set(0); %6
sq.find('Earth Bias 2').set(0); %4
sq.find('Earth Bias 3').set(0); %2
sq.delay(1);

sq.ddsTrigDelay = sq.time;
sq.find('DDS TTL').before(10e-3,1).after(10e-3,0);
sq.find('variable wave plate').set(-2.42);

sq.delay(100e-3);

timeAtDrop = sq.time;

sq.anchor(timeAtDrop);
Traman_pi = 10;
Power_raman_G = 1; %1 %sideband
Power_raman_LG = 1; %0.05 %carrier

% state prep
% sq.dds(1).
Tarp = 10;
t = linspace(0,Tarp,100);
amp = sq.linramp(t,0,1);
sq.dds(1).after(t,110,amp,0);
sq.dds(2).after(t,110,0,0);
% sq.delay(Traman_pi);
sq.dds(1).set(110,0,0);
sq.dds(2).set(110,0,0);

sq.delay(10e-6);
%     sq.find('C6 - N/C').set(1).after(100e-6,0); %SLM trigger

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