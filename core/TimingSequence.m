classdef TimingSequence < handle
    %TimingSequence Defines a sequence of output values
    %based on the events in each channel

    properties
        channels        %Array of TimingControllerChannel objects
        dds             %Array of DDS objects
        ddsTrigDelay    %Offset time between start of sequence and start of DDS
        camDelay        %Time at which to trigger the camera acquisition VI
        
        directory       %Directory where to save sequence builder files
    end

    properties(SetAccess = immutable)
        numDigitalChannels
        numAnalogChannels
        numChannels
    end

    properties(SetAccess = protected)
        data            %Compiled data as a structure with field t (times), d (32-bit unsigned integer array), and a (double-precision 2D array)
        time            %Current time in the object's internal accounting
    end

    properties(Constant)
        SAMPLE_CLK = 1e6;
        MAX_INSTRUCTIONS = 1024;
    end


    methods
        function self = TimingSequence(numDigitalChannels,numAnalogChannels)
            %TimingSequence Constructs a TimingSequence object
            %
            %   sq = TimingSequence(numDigitalChannels,numAnalogChannels)
            %   constructs a TimingSequence object with numDigitalChannels 
            %   digital channels and numAnalogChannels analog channels

            self.numDigitalChannels = numDigitalChannels;
            self.numAnalogChannels = numAnalogChannels;
            self.numChannels = self.numDigitalChannels + self.numAnalogChannels;
            tmp(self.numChannels) = TimingControllerChannel;    %Need to create a temporary variable because self.channels cannot be directly assigned as an array of TimingControllerChannels
            self.channels = tmp;
            for nn = 1:self.numDigitalChannels
                self.channels(nn) = DigitalChannel(nn-1);
            end

            for nn = (self.numDigitalChannels+1):self.numChannels
                self.channels(nn) = AnalogChannel;
            end
            
            self.ddsTrigDelay = 0;
            self.dds = DDSChannel;
            self.dds(2) = DDSChannel;
            for nn = 1:numel(self.dds)
                self.dds(nn).channel = nn;
            end

            self.camDelay = 0;
            
            self.data.t = [];
            self.data.d = [];
            self.data.a = [];
            self.time = 0;
            
            self.directory = 'run-archive';
        end


        function self = reset(self)
            %RESET Resets the TimingSequence channels to their default state
            self.time = 0;
            for nn = 1:self.numChannels
                self.channels(nn).reset;
            end
            for nn=1:numel(self.dds)
                self.dds(nn).reset;
            end
        end

        function ch = digital(self,idx)
            %DIGITAL Returns the digital channel with given index
            %
            %   ch = digital(self) returns an array of all digital
            %   channels
            %   
            %   ch = digital(self,IDX) returns digital channel IDX
            %   This is equivalent to self.channels(IDX);
            
            if nargin < 2
                ch = self.channels(1:self.numDigitalChannels);
            else
                ch = self.channels(idx);
                if ~ch.IS_DIGITAL
                    error('Specified channel %d is not a digital channel',idx);
                end
            end
        end

        function ch = analog(self,idx)
            %ANALOG Returns the analog channel with given index
            %
            %   ch = analog(self) returns an array of all analog channels
            %
            %   ch = analog(self,IDX) returns analog channel IDX
            %   This is equivalent to self.channels(self.numDigitalChannels+IDX);
            
            if nargin < 2
                ch = self.channels((self.numDigitalChannels+1):self.numChannels);
            else
                ch = self.channels(self.numDigitalChannels + idx);
                if ~ch.IS_ANALOG
                    error('Specified channel %d is not a analog channel',idx);
                end
            end
        end

        function ch = find(self,name)
            %FIND Finds a channel with the same name as NAME
            %
            %   ch = self.find(NAME) finds channel ch with name NAME
            ch = [];
            for nn = 1:self.numChannels
                if strcmpi(self.channels(nn).name,name)
                    ch = self.channels(nn);
                    break;
                end
            end
            if isempty(ch)
                for nn = 1:numel(self.dds)
                    if strcmpi(self.dds(nn).name,name)
                        ch = self.dds(nn);
                        break;
                    end
                end
            end
            if isempty(ch)
                error('Channel %s not found.  Check spelling?',name);
            end
        end
        
        function self = anchor(self,time)
            %ANCHOR Sets the latest time for each channel to TIME
            %
            %   sq = anchor(sq,TIME) sets the lastTime property for each channel
            %   to TIME
            
            self.time = time;
            for nn = 1:self.numChannels
                self.channels(nn).anchor(time);
            end
            for nn = 1:numel(self.dds)
                self.dds(nn).anchor(time);
            end
        end

        function self = delay(self,time)
            %DELAY Alias for wait
            self.wait(time);
        end
        
        function self = wait(self,waitTime)
            %WAIT sets latest time for each channel to the current sequence time
            %plus a delay
            %
            %   sq = wait(sq,WAITTIME) sets the lastTime property for each channel
            %   to the current value sq.time + WAITTIME
            self.time = self.time + waitTime;
            self.time = round(self.time*self.SAMPLE_CLK)/self.SAMPLE_CLK;
            for nn = 1:self.numChannels
                self.channels(nn).anchor(self.time);
            end
            for nn = 1:numel(self.dds)
                self.dds(nn).anchor(self.time);
            end
        end

        function self = waitFromLatest(self,waitTime)
            %WAITFROMLATEST sets the lastTime property for each channel to the last
            %chronological time among all channels plus a delay
            %
            %   sq = waitFromLatest(sq,WAITTIME) sets the lastTime property for each
            %   channel to the latest time in the sequence among all channels plus
            %   WAITTIME.
            self.time = self.latest + waitTime;
            for nn = 1:self.numChannels
                self.channels(nn).anchor(self.time);
            end
            for nn = 1:numel(self.dds)
                self.dds(nn).anchor(self.time);
            end
        end

        function time = latest(self)
            %LATEST Returns the latest update time
            %
            %   TIME = latest(sq) returns the latest update time TIME for
            %   sequence sq
            time = 0;
            for nn = 1:self.numChannels
                if self.channels(nn).last >= time
                    time = self.channels(nn).last;
                end
            end
            for nn = 1:numel(self.dds)
                if self.dds(nn).last >= time
                    time = self.dds(nn).last;
                end
            end
            
        end

        function self = check(self)
            %CHECK Applies checks to the sequence to ensure that it doesn't
            %cause the machine to catch fire or otherwise suffer a
            %catastrophic failure
            for nn = 1:numel(self.find('H-bridge Helm').times)
                if self.find('H-Bridge Helm').values(nn) == 1 && self.find('H-Bridge Quad').get(self.find('H-Bridge Helm').times(nn)) == 1
                    error('H-Bridge Helm and H-Bridge Quad CANNOT be on at the same time!');
                end
            end
            ch_ttl = self.find('Redpower TTL');
            ch_pow = self.find('Redpower CW');
%             for nn = 1:numel(ch_ttl.times)
%                 if ch_ttl.values(nn) == 0
%                     ch_pow.at(ch_ttl.times(nn),-0.1);
%                 end
%             end
            
        end

        function r = compile(self)
            %COMPILE Compiles the channel sequences
            %
            %   r = compile(sq) creates a structure that represents the channel data
            %   as seen by the gravimeter control program.  r is equivalent to the internal
            %   property sq.data and has three fields: t, d, and a.  t is an Nx1 array of times
            %   at which updates occur.  d is an Nx1 array of 32-bit unsigned integers that
            %   represents the 32 digital channels.  a is an NxM array of doubles with M the number
            %   of analog channels

            %Check for errors
            self.check;
            %Forms two arrays - one is Nx1 of times, and one is NxNUM_CHANNELS and is values
            t = [];
            v = [];
            for nn = 1:self.numChannels
                self.channels(nn).check.sort;
                [t2,v2] = self.channels(nn).getEvents(true);
                t = [t;t2];   %#ok
                vtmp = NaN(numel(v2),self.numChannels);
                vtmp(:,nn) = v2;
                v = [v;vtmp];   %#ok
            end

            %Now we need to create a smaller list of times and values which correspond
            %to unique updates
            [t,k] = sort(round(t*self.SAMPLE_CLK));
            v = v(k,:);
            buf = NaN(size(t,1),1+self.numChannels);    %One column for time, one for all other channels
            buf(1,:) = [t(1),v(1,:)];
            numBuf = 1;
            for nn = 2:numel(t)
                if t(nn) ~= t(nn-1)
                    numBuf = numBuf + 1;
                    buf(numBuf,1) = t(nn);
                end
                idx = find(~isnan(v(nn,:)));
                buf(numBuf,1+idx) = v(nn,idx);
            end
            buf = buf(1:numBuf,:);
            
            %Replace NaNs with previous values
            for nn = 2:size(buf,1)
                tmpOld = buf(nn-1,2:end);
                tmpCurrent = buf(nn,2:end);
                tmpCurrent(isnan(tmpCurrent)) = tmpOld(isnan(tmpCurrent));
                buf(nn,2:end) = tmpCurrent;
                
                %Last two digital channels switch at each update
                buf(nn,1+self.numDigitalChannels+(-1:0)) = double(~buf(nn-1,1+self.numDigitalChannels+(-1:0)));
            end
            %One column for time, one for all digital channels, and one for each analog channel
            bits = zeros(1,self.numDigitalChannels);
            for nn = 1:numel(bits)
                bits(nn) = self.channels(nn).bit;
            end
            
            if size(buf,1) > self.MAX_INSTRUCTIONS
                error('Instruction set size of %d is larger than maximum size of %d!',size(buf,1),self.MAX_INSTRUCTIONS);
            end
            
            self.data.t = buf(:,1)/self.SAMPLE_CLK;
            self.data.d = uint32(sum(buf(:,1+(1:self.numDigitalChannels)).*repmat(2.^bits,size(buf,1),1),2));
            self.data.a = buf(:,1+((self.numDigitalChannels+1):self.numChannels));
            
            for nn = 1:numel(self.dds)
                self.data.dds(nn) = self.dds(nn).compile(self.ddsTrigDelay);
            end

            self.data.camDelay = self.camDelay;
            
            r = self.data;
        end

        function self = reduce(self)
            %REDUCE reduces updates to only those which change channel values

            for ch = self.channels
                ch.reduce;
            end
        end

        function self = loadCompiledData(self,data)
            %LOADCOMPILEDDATA Loads compiled data into channel values
            %
            %   sq = loadCompiledData(sq,DATA) loads compiled data DATA
            %   into the sequence structure.  DATA should have fields t, d,
            %   and a.  t should be Nx1 double, d should be Nx1 uint32, and 
            %   a should be NxM double with M the number of analog channels

            self.reset;
            for nn = 1:numel(data.t)
                for ch = self.digital()
                    ch.at(data.t(nn),bitget(data.d(nn),ch.bit+1));
                end

                for mm = 1:self.numAnalogChannels
                    self.analog(mm).at(data.t(nn),data.a(nn,mm));
                end
            end
        end

        function plot(self,offset)
            %PLOT Plots all channel sequences
            %
            %   tc.plot plots all channel sequences on the same graph
            %
            %   tc.plot(offset) plots all channel sequences on the
            %   same graph but with each channel's sequence offset from the
            %   next by offset
            jj = 1;
            if nargin < 2
                offset = 0;
            end
            str = {};
            for nn = 1:self.numChannels
                if self.channels(nn).exists
                    self.channels(nn).plot('offset',(jj-1)*offset,'finaltime',self.latest);
                    hold on;
                    if isempty(self.channels(nn).name)
                        str{jj} = sprintf('Ch %d',nn);  %#ok<AGROW>
                    else
                        str{jj} = sprintf('%s',self.channels(nn).name); %#ok<AGROW>
                    end
                    jj = jj+1;
                end
            end
            hold off;
            legend(str);
            xlabel('Time [s]');
        end
        
        function disp(self)
            %DISP Creates a display of the object
            fprintf(1,'\tTimingSequence object with properties:\n');
            fprintf(1,'\t\tNumber of Digital Channels: %d\n',self.numDigitalChannels);
            fprintf(1,'\t\t Number of Analog Channels: %d\n',self.numAnalogChannels);
            fprintf(1,'\t\t      Instruction set size: %d\n',size(self.data.t,1));
            fprintf(1,'\n');
            fprintf(1,'\tDigital channels:\n');
            fprintf(1,'\t%-10s\t\t%-20s\t\t%-s\n','Port','Name','Description');
            for ch = self.digital()
                if ~isempty(ch.name)
                    fprintf(1,'\t%-10s\t\t%-20s\t\t%-s\n',ch.port,ch.name,ch.description);
                end
            end
            fprintf(1,'\n');
            fprintf(1,'\tAnalog channels:\n');
            fprintf(1,'\t%-10s\t\t%-20s\t\t%-s\n','Port','Name','Description');
            for ch = self.analog()
                if ~isempty(ch.name)
                    fprintf(1,'\t%-10s\t\t%-20s\t\t%-s\n',ch.port,ch.name,ch.description);
                end
            end
        end

    end

    methods(Static)
        function v = minjerk(t,vi,vf)
            T = max(t)-min(t);
            t = (t-min(t))./T;
            v = vi + (vf-vi).*(10*t.^3-15*t.^4+6*t.^5);
        end
        
        function v = linramp(t,vi,vf)
            v = vi+(vf-vi)/(t(end)-t(1))*(t-t(1));
        end
        
        function v = expramp(t,vi,vf,T)
            e = exp(-(t-t(1))/T);
            v = vi*e+(vf-vi*e(end))./(1-e(end)).*(1-e);
        end

        function sq = buildFromCompiledData(data)
            sq = TimingSequence(32,size(data.a,2));
            sq.loadCompiledData(data);
        end
    end


end