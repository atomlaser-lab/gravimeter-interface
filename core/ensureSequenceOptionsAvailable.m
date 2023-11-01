function varargout = ensureSequenceOptionsAvailable(varargin)
    % Ensure SequenceOptions object is available in the base workspace
    % ...

    % Get all variable names in the workspace
    allVarNames = evalin('base', 'who');

    % Check if any variable is an instance of SequenceOptions
    seqOptVarIndices = cellfun(@(varName) isa(evalin('base', varName), 'SequenceOptions'), allVarNames);
    seqOptVarNames = allVarNames(seqOptVarIndices);

    if isempty(seqOptVarNames)
        if nargin == 1 && ~exist(varargin{1},'var')
            error('define the sequence option variable name first');
        end
        % Create a new instance of SequenceOptions if it does not exist
        opt = SequenceOptions;
        cprintf('Keywords','Initialising Sequence Options\n')
        assignin('base', 'opt', opt);

    else
        if numel(seqOptVarNames) > 1
            error('Multiple SequenceOptions found in the workspace');
        end
        % Use the existing instance of SequenceOptions
        opt = evalin('base', seqOptVarNames{1});
    end

    % Check input arguments and update options accordingly
    if nargin == 1
        if ~isa(varargin{1}, 'SequenceOptions')
            error('If using only one argument, it must be of type SequenceOptions');
        end
    elseif nargin > 1
        if ~isa(varargin{1}, 'SequenceOptions')
            error('First argument must be of type SequenceOptions');
        end
        opt.replace(varargin{1});
        opt.set(varargin{2:end});
    end

    if nargout == 0
        % Search for the opt variable in the base workspace
        varExists = evalin('base', 'exist(''opt'', ''var'')');
        if varExists
            opt = evalin('base', 'opt');
        end
    else
        varargout{1} = opt;
    end
end
