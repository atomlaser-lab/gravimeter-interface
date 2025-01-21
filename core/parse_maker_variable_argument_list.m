function opt = parse_maker_variable_argument_list(varargin)

opt = SequenceOptions('load_time',15,'detuning',0,'tof',20e-3,'redpower',2,...
    'raycus',2);

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
else
    error('Either supply a single SequenceOptions argument, or supply a set of name/value pairs, or supply a SequenceOptions argument followed by name/value pairs');
end