classdef RolloverCounter < handle
    %ROLLOVERCOUNTER Class defining multiple counters which increment
    %sequentially - just like a car odometer!
    properties
        N           %Number of indices
        i           %Array of current indices
        
        initial     %Inital value of each index
        final       %Final value of each index
    end
    
    methods
        function self = RolloverCounter(varargin)
            %ROLLOVERCOUNTER Creates a RolloverCounter object
            %
            %   counter = RolloverCounter(FINAL) creates an object with
            %   final index values specified by FINAL.  The number of
            %   elements in the vector FINAL specifies the number of
            %   indices to count, with FINAL(1) being the final index value
            %   for the first counter, FINAL(2) the second counter, and so
            %   on.
            %
            %   counter = RolloverCounter(FINAL,INITIAL) creates an object
            %   with final index values given by FINAL and initial values
            %   for each index given by INITIAL.  FINAL and INITIAL must
            %   have the same number of elements.
            if nargin > 0
                self.setup(varargin{:});
            else
                self.initial = 1;
                self.final = Inf;
                self.N = 1;
            end
        end
        
        function self = setup(self,varargin)
            %SETUP Sets up the RolloverCounter object COUNTER
            %
            %   COUNTER = COUNTER.setup(FINAL) sets  the final index values
            %   for COUNTER as specified by FINAL.  The number of elements
            %   in the vector FINAL specifies the number of indices to
            %   count, with FINAL(1) being the final index value for the
            %   first counter, FINAL(2) the second counter, and so on.
            %
            %   COUNTER = COUNTER.setup(FINAL,INITIAL) creates an object
            %   with final index values given by FINAL and initial values
            %   for each index given by INITIAL.  FINAL and INITIAL must
            %   have the same number of elements.
            %
            %   COUNTER = COUNTER.setup('var',v1,v2,...) creates an object
            %   with N set to the number of variable arguments v1,v2,...
            %   and final values set to the number of elements in each
            %   variable argument v1,v2,...
            if numel(varargin) == 1
                final = varargin{1}; %#ok<*PROPLC>
                initial = ones(size(final));
            elseif ~ischar(varargin{1})
                final = varargin{1};
                initial = varargin{2};
            else
                final = ones(numel(varargin(2:end)),1);
                for nn = 1:numel(varargin(2:end))
                    final(nn) = numel(varargin{nn+1});
                end
                initial = ones(size(final));
            end
            
            if numel(initial) ~= numel(final)
                error('Initial and final indices must have the same number of elements!');
            end
            
            self.N = numel(initial);
            self.initial = initial(:)';
            self.final = final(:)';
            self.reset;
        end
        
        function self = reset(self)
            %RESET Resets the counter object COUNTER to inital values
            %
            %   COUNTER = COUNTER.reset() Resets the COUNTER object
            self.i = self.initial;
        end
        
        function im = imax(self,idx)
            %IMAX Returns the number of runs for each index
            %
            %   IM = COUNTER.IMAX()
            %
            if nargin == 1
                idx = 1;
            end
            im = self.final(idx) - self.initial(idx) + 1;
        end
        
        function r = total(self)
            %TOTAL Returns the total number of runs
            %
            %   R = COUNTER.total() returns the total number of runs for
            %   COUNTER in R
            r = prod(self.imax(1:self.N));
        end
        
        function c = current(self)
            %CURRENT Returns the current, flattened index
            %
            %   C = COUNTER.current() returns the current, flattened index
            %   in C for COUNTER.  This value of C/COUNTER.total() is the
            %   completion rate of the counter
            c = 1;
            im = self.imax(1:self.N);
            for nn = 1:self.N
                c = c + (self.i(nn)-1)*prod(im((nn-1):-1:1));
            end
        end
        
        function c = now(self)
            %NOW Alias of CURRENT()
            c = self.current;
        end
        
        function self = increment(self)
            %INCREMENT Increments the counter by 1
            %
            %   COUNTER = COUNTER.increment() increments the counter by 1,
            %   handling roll-over as necessary
            for nn = 1:self.N
                if nn == 1
                    self.i(nn) = self.i(nn) + 1;
                end
                
                if self.i(nn) > self.final(nn)
                    self.i(nn) = self.initial(nn);
                    if nn < self.N
                        self.i(nn+1) = self.i(nn+1) + 1;
                    end
                end
                    
            end
        end

        function self = decrement(self)
            %DECREMENT Decrements the counter by 1
            %
            %   COUNTER = COUNTER.DECREMENT() decrements the counter by 1,
            %   handling roll-over as necessary
            for nn = 1:self.N
                if nn == 1
                    self.i(nn) = self.i(nn) - 1;
                end

                if self.i(nn) < self.initial(nn)
                    self.i(nn) = self.final(nn);
                    if nn < self.N
                        self.i(nn+1) = self.i(nn+1) - 1;
                    end
                end
            end
        end
        
        function r = done(self,idx)
            %DONE Indicates if counter is done
            %
            %   R = COUNTER.done() returns true if COUNTER is on its last
            %   value, false otherwise
            %
            %   R = COUNTER.done(IDX) returns true if the counter for index
            %   IDX is on its last value
            if nargin < 2
                r = (self.current() == self.total());
            else
                r = self.i(idx) == self.final(idx);
            end
        end
        
        function varargout = print(self)
            %PRINT Prints a summary of the current state of the counter
            %
            %   COUNTER.print() prints the summary to the command line
            %
            %   S = COUNTER.print() creates a string S with the summary
            s = '';
            for nn = 1:self.N
                s = [s,sprintf('%d/%d',self.i(nn),self.imax(nn))]; %#ok<*AGROW>
                if nn ~= self.N
                    s = [s,' '];
                end
            end
            if nargout == 0
                fprintf(1,'%s\n',s);
            else
                varargout{1} = s;
            end
        end
        
        function B = subsref(self,S)
            switch S(1).type
                case '.'
                    B = builtin('subsref',self,S);
                case '()'
                    if length(S) < 2
                        B = builtin('subsref',self.i,S);
                    else
                        B = builtin('subsref',self,S);
                    end
                case '{}'
                    error('Not a supported reference');
            end
        end
    end
    
    
end