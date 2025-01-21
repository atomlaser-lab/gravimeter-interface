classdef SequenceOptionsAbstract < matlab.mixin.Copyable
    %SEQUENCEOPTIONSABSTRACT Defines an abstract class for passing options
    %to the gravimeter make sequence function
    
    
    methods(Abstract)
        self = setDefaults(self)
    end

    methods       
        function self = set(self,varargin)
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            else
                for nn = 1:2:numel(varargin)
                    try
                        if isa(self.(varargin{nn}),'SequenceOptionsAbstract')
                            self.(varargin{nn}).set(varargin{nn+1}{:});
                        else
                            self.(varargin{nn}) = varargin{nn+1};
                        end
                    catch
                    end
                end
            end
        end

        function self = replace(self,opt)
            p = properties(opt);
            for nn = 1:numel(p)
                if ~isempty(opt.(p{nn}))
                    if isa(self.(p{nn}),'SequenceOptionsAbstract')
                        self.(p{nn}).replace(opt.(p{nn}));
                    else
                        self.(p{nn}) = opt.(p{nn});
                    end
                end
            end
        end

        function s = print(self)
            p = properties(self);
            sargs = {};
            for nn = 1:numel(p)
                v = self.(p{nn});
                if isempty(v)
                    sargs{nn} = sprintf('''%s'',%s',p{nn},'[]'); %#ok<*AGROW> 
                elseif ischar(v) || isstring(v)
                    sargs{nn} = sprintf('''%s'',''%s''',p{nn},v);
                elseif isa(v,'SequenceOptionsAbstract')
                    sargs{nn} = v.print;
                    sargs{nn}(end) = '';
                    sargs{nn} = regexprep(sargs{nn},'\n','\n    ');
                    sargs{nn} = sprintf('''%s'',%s',p{nn},sargs{nn});
                else
                    sargs{nn} = sprintf('''%s'',%.6g',p{nn},v);
                end
            end
            s = strjoin(sargs,',...\n    ');
            s = sprintf('%s(%s);',class(self),s);
        end
        
    end
    
end