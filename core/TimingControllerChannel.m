classdef TimingControllerChannel < handle & matlab.mixin.Heterogeneous
    %TimingControllerChannel Defines a generic timing
    %controller channel.  This class should be extended
    %before it is used.
    %
    %   The matlab.mixin.Heterogeneous is needed to allow for arrays
    %   of different subclasses of this base class
    properties
        name                        %Name of the channel
        port                        %Fixed port of the channel
        description                 %Description of the channel
        conversion_function         %Function handle converting "real" values to volts
        units                       %The "real" units for this channel
        manual                      %Manual value
    end
    
    properties(SetAccess = protected)
        default     %Default value for this channel
%         numSubValues%Number of sub values, expressed as columns in values property
        
        values      %Array of values in channel sequence.
        times       %Array of times in the channel sequence in seconds
        numValues   %Number of time/value pairs
        
        lastTime    %Last time written - used for before/after functions
        
        bounds      %Allowed bounds for values
        IS_DIGITAL  %Indicates if a channel is a digital channel
        IS_ANALOG   %Indicates if a channel is an analog channel
        IS_DDS      %Indicates if a channel is a DDS channel
    end
    
    methods
        function ch = TimingControllerChannel
            %TimingControllerChannel Contructs a channel

            ch.IS_DIGITAL = false;
            ch.IS_ANALOG = false;
            ch.IS_DDS = false;
            ch.bounds = [0,0];
            ch.default = 0;
            ch.manual = ch.default;
            ch.numValues = 0;
            ch.lastTime = 0;
            ch.port = '';
            ch.description = '';
            ch.times = [];
            ch.values = [];
            ch.setConversionFunction(@(x) x,'V');
        end
        
        function ch = setName(ch,name,port,description)
            %SETNAME Sets the name and optionally the description
            %
            %   ch = setName(ch,NAME,PORT) sets the name property to NAME and
            %   port property to PORT
            %
            %   ch = setName(ch,NAME,PORT,DESC) sets the name property to NAME
            %   and the description property to DESC
            ch.name = name;
            if nargin >= 3
                ch.port = port;
            end
            if nargin >= 4
                ch.description = description;
            end
        end

        function ch = setDefault(ch,v)
            %SETDEFAULT Sets the default value
            %
            %   ch = setDefault(ch,v) sets the default value of channel ch
            %   to v
            ch.checkValue(v);
            ch.default = v;
            if numel(ch.times) == 0
                ch.times = 0;
                ch.values = ch.default;
            else
                idx = find(ch.times == 0);
                if isempty(idx)
                    ch.times(end+1) = 0;
                    ch.values(end+1,:) = ch.default;
                    lt = ch.lastTime;
                    ch.sort;
                    ch.lastTime = lt;
                else
                    ch.values(idx,:) = ch.default;
                end
            end
            ch.numValues = numel(ch.times);
        end
        
        function [t,v] = getEvents(ch,get_volts)
            %getEvents Returns the times and values as separate Nx1 arrays.
            % 
            %   Events are checked for errors and sorted before being
            %   returned.
            %   [t,v] = ch.getEvents returns times t and values v
            %   Events are always returned starting at time 0 - if no 
            %   value at time 0 is specified, the default value is used
            if nargin < 2
                get_volts = 0;
            end
            ch.check;
            ch.sort;
            if isempty(ch.numValues) || (ch.numValues == 0)
                t = 0;
                v = ch.default;
            elseif ch.times(1) == 0
                t = ch.times;
                v = ch.values;
            else
                t = [0;ch.times];
                v = [ch.default;ch.values];
            end
            if get_volts
                v = ch.convert(v);
            end
        end
        
        function r = exists(ch)
            %EXISTS Indicates if channel sequence is more than just the
            %default value
            
            r = ~( (numel(ch.times) == 1 && ch.times(1) == 0) || isempty(ch.times) );

        end

        function ch = setBounds(ch,bounds)
            %SETBOUNDS Sets the bounds for channel
            %
            %   CH = setBounds(CH,BOUNDS) Sets the bounds to the 2-element
            %   array BOUNDS.  BOUNDS can be in any order
            ch.bounds(1) = min(bounds);
            ch.bounds(2) = max(bounds);
        end

        function ch = setConversionFunction(ch,func_in,units)
            %SETCONVERSIONFUNCTION Sets the conversion function from "real"
            %values (MHz, G/cm, etc) to volts
            %
            %   CH = setConversionFunction(CH,FUNC_IN,UNITS) sets
            %   conversion function to FUNC_IN associated with conversion
            %   from UNITS to Volts
            ch.conversion_function = func_in;
            if nargin > 2
                ch.units = units;
            end
        end

        function values_out = convert(ch,values_in)
            %CONVERT Converts "real" values (like MHz, G/cm, etc) to volts
            %
            %   VALUES_OUT = convert(CH,VALUES_IN) converts "reaL'VALUES_IN
            %   to VALUES_OUT in volts using CH.CONVERSION_FUNCTION
            values_out = ch.conversion_function(values_in);
        end
        
        function ch = at(ch,time,value,varargin)
            %AT Adds a value at the given time
            %
            %   ch = ch.at(TIME,VALUE) adds VALUE to the events at the time 
            %   given by TIME.  Assumes that TIME and VALUE are scalars.  TIME must 
            %   be in seconds.  Sets the lastTime property to TIME
            %
            %   CH = CH.at(TIME,VALUE) if TIME and VALUE are Nx1 arrays
            %   adds each row in VALUE to events at time given by
            %   corresponding row in VALUE. This uses a recursive call to
            %   AT, so may run into memory issues 
            %
            %   CH = CH.at(TIME,VALUE) if TIME is an Nx1 array and VALUE is
            %   an NxM array adds each row in VALUE to events at time given
            %   by corresponding row in VALUE. This uses a recursive call
            %   to AT, so may run into memory issues
            %
            %   CH = CH.at(TIME,VALUEFUNC) if TIME is an Nx1 array and
            %   VALUEFUNC is a function handle, adds the calculated values
            %   VALUEFUNC(TIME) to events.  Uses recursive calls, so may
            %   run into memory issues
            %
            %   CH = CH.at(TIME,VARARGIN) if time is an Nx1 array and each
            %   of the VARARGIN elements are Nx1 array, adds the values for
            %   each element to the corresponding row in VALUE. This uses a
            %   recursive call to AT, so may run into memory issues

            time = time(:);
            N = numel(time);
            if numel(ch) > 1
                %If an array of channels is passed, loop through each individually
                for nn = 1:numel(ch)
                    ch(nn).at(time,value,varargin{:});
                end
            elseif ~isa(value,'function_handle') && (N > 1) && (numel(varargin) ~= 0)
                %If TIME is Nx1 and VALUE and VARARGIN{...} are each
                %Nx1, recursively add events
                value = value(:);
                Nv = numel(value);
                if Nv == 1
                    value = value*ones(N,1);
                end
                for nn = 1:numel(varargin)
                    tmp = varargin{nn};
                    value(:,nn+1) = tmp(:).*ones(N,1);
                end
                ch.at(time,value);
            elseif isa(value,'function_handle') && (N > 1)
                %If TIME is an array and VALUE is a function handle, loop through each time and calculate a value
                v = value(time);
                ch.at(time,v);
            else
                %Otherwise add single events
                if numel(varargin) ~= 0
                    value = value(:);
                    for nn = 1:numel(varargin)
                        value(:,nn+1) = varargin{nn};
                    end
                elseif any(size(value) == 1)
                    value = value(:);
                end
                ch.checkValue(value);   %Check that value is within bounds

                time = round(time*TimingSequence.SAMPLE_CLK)/TimingSequence.SAMPLE_CLK;   %Round time to multiple of sample clock
                for nn = 1:numel(time)
                    idx = (ch.times == time(nn));
                    if sum(idx) == 0
                        %If this has not been used previously, add a new
                        %value at this time
                        N = ch.numValues + 1;
                        ch.values(N,:) = value(nn,:);
                        ch.times(N,1) = time(nn);
                        ch.numValues = N;
                        ch.lastTime = time(nn);
                    else
                        %If time has been used, replaced that value
                        ch.values(idx,:) = value(nn,:);
                        ch.times(idx) = time(nn);
                        ch.lastTime = time(nn);
                    end
                end
%                 idx = find(ch.times==time,1,'first');   %Try and find first time
%                 if isempty(idx)
%                     %If this time has not been used previously, add a new value at this time
%                     N = ch.numValues+1;
%                     ch.values(N,:) = value;
%                     ch.times(N,1) = time;
%                     ch.numValues = N;
%                     ch.lastTime = time;
%                 else
%                     %If time has been used, replace that value
%                     ch.values(idx,:) = value;
%                     ch.times(idx,1) = time;
%                     ch.lastTime = time;
%                 end
            end
        end
        
        function ch = on(ch,varargin)
            %ON Alias of AT method
            ch.at(varargin{:});
        end
        
        function ch = set(ch,value,varargin)
            %SET Sets a value at the current lastTime
            %
            %   ch = set(ch,value) sets the current value to be value at time
            %   lastTime
            if numel(ch) > 1
                for nn = 1:numel(ch)
                    ch(nn).at(ch(nn).lastTime,value,varargin{:});
                end
            else
                ch.at(ch.lastTime,value,varargin{:});
            end
        end
        
        function ch = after(ch,delay,value,varargin)
            %AFTER Adds a value to the events after the last added event
            %
            %   ch = ch.after(DELAY,VALUE) adds value to the events a time
            %   DELAY seconds after the property lastTime.  Note that this
            %   is not necessarily the latest time in the sequence.
            
            if numel(ch) > 1
                for nn = 1:numel(ch)
                    ch(nn).after(delay,value,varargin{:});
                end
            else
                time = ch.lastTime+delay;
                ch.at(time,value,varargin{:});
            end
        end
        
        function ch = before(ch,delay,value,varargin)
            %BEFORE Adds a value to the events before the last added event
            %
            %   ch = ch.before(DELAY,VALUE) adds value to the events a time
            %   DELAY seconds before the property lastTime.  Note that this
            %   is not necessarily the latest time in the sequence.
            
            if numel(ch) > 1
                for nn = 1:numel(ch)
                    ch(nn).before(delay,value,varargin{:});
                end
            else
                time = ch.lastTime-delay;
                ch.at(time,value,varargin{:});
            end
        end
        
        function ch = anchor(ch,time)
            %ANCHOR Sets the lastTime property
            %
            %   ch.anchor(TIME) sets the lastTime property to TIME
            
            if numel(ch) > 1
                for nn = 1:numel(ch)
                    ch(nn).anchor(time);
                end
            else
                ch.lastTime = round(time*TimingSequence.SAMPLE_CLK)/TimingSequence.SAMPLE_CLK;
            end
        end
        
        function [time,value] = last(ch)
            %LAST Returns the last time and last value
            %
            %   [t,v] = ch.last returns the last time t and last value v
            [t,v] = ch.getEvents;
            time = t(end);
            if any(size(v) ~= 1)
                value = v(end,:);
            else
                value = v(end);
            end
        end
        
        function ch = reset(ch)
            %RESET Resets the channel sequence so that there are no events
            %   ch = ch.reset resets the channel
            ch.times = [];
            ch.values = [];
            ch.setDefault(ch.default);
            ch.sort;
        end
        
        function ch = sort(ch)
            %SORT Sorts the events so that they are ordered chronologically
            %
            %   ch = ch.sort sorts the events.  The lastTime property is
            %   set to the last time in the sorted events
            if numel(ch.times)>0
                [B,K] = sort(ch.times);
                ch.times = B;
                ch.values = ch.values(K,:);
                ch.lastTime = ch.times(end);
            end
            ch.numValues = numel(ch.times);
        end
        
        function ch = checkTimes(ch)
            %CHECKTIMES Checks times to make sure that they are all >= 0
            %
            %   Also checks to see if unique events actually occur, and
            %   removes sequence if nothing happens
            %
            %   ch = ch.check checks the event times and removes sequence
            %   if nothing happens
%             if numel(unique(ch.values))==0
%                 ch.reset;
%             end
            if any(ch.times<0)
                error('All times must be greater than 0 (no acausal events)!');
            end
        end

        function ch = checkValue(ch,v)
            %CHECKVALUE Checks a given value to make sure it is in a valid range
            %
            %   CH = CHECKVALUE(CH,V) checks value v to make sure it is within 
            %   a valid range
            if any(v < ch.bounds(1)) || any(v > ch.bounds(2))
                error('Value outside of bounds [%.3f,%.3f]',ch.bounds(1),ch.bounds(2));
            end
        end

        function ch = check(ch)
            %CHECK Checks events to ensure correct times and values are given
            %
            %   CH = CHECK(CH) Checks events to make sure times are positive and that unique times
            %   are given.  Also checks that values are within bounds
            ch.checkTimes();
            for nn = 1:ch.numValues
                ch.checkValue(ch.values(nn,:));
            end
        end

        function ch = reduce(ch)
            %REDUCE Reduces channel events to only those which change the value

            ch.sort;
            [t,v] = ch.getEvents;
            ch.times = t(1);ch.values = v(1);
            idx = 2;
            for nn = 2:numel(v)
                if v(nn) ~= v(nn-1)
                    ch.times(idx,1) = t(nn);
                    ch.values(idx,:) = v(nn,:);
                    idx = idx + 1;
                end
            end
            ch.sort;
        end

        function v = get(self,t)
            %GET Returns the channel value at a given time
            %
            %   V = GET(T) returns the channel value at a time T
            self.sort;
            idx = find(self.times < t,1,"last");
            v = self.values(idx);
        end
        
        function print(ch,idx)
            if nargin < 2
                idx = 1:numel(ch.times);
            end
            for nn = idx
                fprintf(1,['Time: % 12.6f, Value: ',repmat('% 12.3f ',1,size(ch.values,2)),'%s','\n'],ch.times(nn),ch.values(nn,:),ch.units);
            end 
        end

        
        function varargout = plot(ch,varargin)
            %PLOT Plots the current sequence as a function of time.
            %
            %   ch.plot plots the current sequence as a function of time.
            %   If there are no events, a message is displayed.
            %
            %   ch.plot(AX) plots the current sequence on the axes given by
            %   AX.
            %
            %   ch.plot(AX,'Name','Value',...) plots the current sequnce on
            %   axes AX with properties given by Name and Value.  Name can
            %   be OFFSET, which gives a vertical offset on the channel,
            %   FINALTIME, which plots the sequence out to FINALTIME, and
            %   PLOTARGS, which is a cell array of plot-arguments
            %
            %   ch.plot('Name','Value') plots as above in the current axes
            %
            if numel(varargin) >= 1 && all(ishandle(varargin{1})) && strcmpi(get(varargin{1},'type'),'axes')
                ax = varargin{1};
                varargin = varargin(2:end);
            else
                ax = gca;
            end
            
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            else
                returnHandle = false;
                plotargs = {};
                plotIdx = 1:size(ch.values,2);
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case 'plotargs'
                            plotargs = v;
                        case 'returnhandle'
                            returnHandle = v;
                        case 'plotidx'
                            plotIdx = v;
                    end
                end
            end
            
            if ~ch.exists && ~returnHandle
                if nargout > 0
                    varargout{1} = [];
                end
                return
            end
            
            [tplot,vplot] = ch.getPlotValues(varargin{:});
            cargs = {'linewidth',1.5,'tag',ch.name};
            if numel(plotargs) > 0
                cargs = [cargs,plotargs];
            end
            if ~ch.IS_DIGITAL
                h = plot(ax,tplot,vplot,'.-',cargs{:});
            else
                h = plot(ax,tplot,vplot,'.--',cargs{:});
            end
            
            if nargout > 0
                varargout{1} = h;
            end
        end
        
        function [tplot,vplot] = getPlotValues(ch,varargin)
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be in name/value pairs');
            else
                offset = 0;
                finalTime = [];
                returnHandle = false;
                plotIdx = 1:size(ch.values,2);
                plot_volts = 0;
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case 'offset'
                            offset = v;
                        case 'finaltime'
                            finalTime = v;
                        case 'returnhandle'
                            returnHandle = v;
                        case 'plotidx'
                            plotIdx = v;
                        case 'plotvolts'
                            plot_volts = 1;
                    end
                end
            end
            [t,v] = ch.getEvents(plot_volts);
            if ~ch.exists && ~returnHandle
                tplot = [];
                vplot = [];
                return
            end
            
            if t(end) ~= finalTime
                t = [t;finalTime];
                v = [v;v(end,:)];
            end
            tplot = sort([t;t-1/TimingSequence.SAMPLE_CLK]);
            tplot = tplot(tplot >= 0);
            vplot = zeros(numel(tplot),numel(plotIdx));
            for nn = 1:numel(plotIdx)
                vplot(:,nn) = interp1(t,v(:,plotIdx(nn)),tplot,'previous');
            end
            vplot = vplot + offset;
        end
        
    end
    
    
end