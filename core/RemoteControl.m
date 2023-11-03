classdef RemoteControl < handle
    properties
        %% TCPIP properties
        conn            %TCPIP connection
        connected       %Is LabVIEW client connected?
        %% Sequence properties
        status          %Current status of run: RUNNING or STOPPED
        sq              %Sequence object representing current sequence
        makerCallback   %Callback function for creating a TimingSequence object
        %% Data properties
        mode            %Mode of callback function: SET, ANALYZE, or INIT
        devices         %Structure listing MATLAB devices used in callback
        data            %Data structure to use in callback function
        callback        %Callback function, takes argument of Rebeka object
        %% DDS properties
        mog             %MOGLabs parent object
    end
    
    properties(SetAccess = protected)
        remoteAddress = 'localhost';  %Connect to local host
        remotePort = 6666;            %Remote port to use
    end %end constant properties
    
    properties(SetAccess = immutable)
        c               %Rollover counter object, keeps track of runs
        err             %RemoteControlErrorHandler object, used for handling errors
    end
    
    properties(Constant, Hidden=true)
        readyWord = 'ready';          %Word indicating that client is ready
        startWord = 'start';          %Word telling host to start
        endWord = 'end';              %Word telling host to stop TCP loop
        uploadDWord = 'uploadD';      %Word telling host to upload digital (uint32) data
        uploadAWord = 'uploadA';      %Word telling host to upload analog (float) data
        uploadCamDelay = 'cam-delay'; %Word telling host to set the camera acquisition delay

        SET = 'set/check';
        ANALYZE = 'analyze';
        INIT = 'init';
        
        RUNNING = 'running';
        STOPPED = 'stopped';
    end    
    
    events
        sequenceChanged
    end
    
    methods
        function self = RemoteControl(varargin)
            self.setRemoteProperties(varargin{:});
            self.connected = false;
            self.mode = self.INIT;
            self.status = self.STOPPED;
            self.makerCallback = @makeBEC_scratch;
            self.c = RolloverCounter();
            self.err = RemoteControlErrorHandler;
            self.reset;
        end %end constructor
        
        function self = setRemoteProperties(self,varargin)
            if numel(varargin) >= 1
                self.remoteAddress = varargin{1};
            end
            if numel(varargin) >= 2
                self.remotePort = varargin{2};
            end
        end
        
        function open(self)
            %OPEN Opens a tcpip port           
            %open Creates and opens a TCP conn.  Waits for ready word
            if isempty(self.conn)
                fprintf(1,'Attempting connection...\n');
                self.conn = tcpclient(self.remoteAddress,self.remotePort);
                self.conn.configureTerminator('CR/LF');
                self.conn.configureCallback('terminator',@(src,event) self.resp(src,event))
                self.conn.OutputBufferSize = 2^24;
                self.conn.InputBufferSize = 2^24;
                fprintf(1,'Connection successful!\n');
                self.connected = true;
            end
        end %end open
        
        function setFunc(self,cb)
            %SETFUNC Sets the BytesAvailableFcn to self.resp()
            self.open;
            if nargin < 2
                self.conn.configureCallback('terminator',@(src,event) self.resp(src,event));
            else
                self.conn.configureCallback('terminator',@(~,~) cb());
            end
        end
        
        function r = read(self)
            %READ Reads available data from TCP connection
            r = self.conn.readline();
        end %end read
        
        function stop(self)
            %STOP Releases client from remote control and closes TCP conn
            self.conn = [];
            fprintf(1,'Remote control session terminated\n');
            self.connected = false;
            self.status = self.STOPPED;
        end
        
        function delete(self)
            %DELETE Deletes this object
            %
            %   Closes then deletes the tcpip connection with the LabVIEW
            %   interface before deleting the object
            self.stop;
        end

        function self = make(self,varargin)
            %MAKE Makes the sequence to be uploaded
            %
            %   r = make(r,varargin) runs r.sq = r.makerCallback(varargin{:})
            %   If r.makerCallback is empty, then uses default makeSequence() function
            if isempty(self.makerCallback) || ~isa(self.makerCallback,'function_handle')
                self.makerCallback = @makeSequence;
            end
            self.sq = self.makerCallback(varargin{:});
            notify(self,'sequenceChanged');
        end
        
        function self = upload(self,data)
            %UPLOAD uploads data to host
            %
            %   r = upload(r) uploads data to control interface using the 
            %   current sequence stored in the r.sq field
            %
            %   r = upload(r,data) with r the RemoteControl object and data a
            %   2D array with times in the first column, a 32 bit digital
            %   value in the second column, and 24 analog values in the rest
            if nargin < 2
                data = self.sq.compile;
            end

            if isnumeric(data)
                if size(data,2) ~= 26
                    error('Numeric input array must have 26 columns!');
                end
                d = uint32(round(data(:,2)));
                a = data(:,[1,3:end]);
            elseif isstruct(data)
                d = uint32(data.d);
                a = [data.t,data.a];
                if size(d,1) ~= size(a,1)
                    error('Analog and digital columns must have the same size!');
                elseif size(data.a) ~= 24
                    error('Data ''a'' field must have 24 columns');
                end
            end
            
            %% Upload DDS data
            if strcmpi(data.dds(1).power_conversion_method,DDSChannel.POWER_CONVERSION_HEX_INTERP)
                self.uploadBinaryDDSData(data.dds);
            else
                self.uploadDDSData(data.dds);
            end
            %% Upload R&S synthesizer list data if present
            if isfield(self.devices,'rs') && isfield(self.devices.rs,'list') && ~isempty(self.devices.rs.list.freq)
                self.devices.rs.writeList;
            end
            
            %% Open connection with LabVIEW VI and set camera acquisition delay
            self.open;
            self.conn.writeline(self.uploadCamDelay);
            s = sprintf('%.1f',data.camDelay);
            pause(0.1);
            self.conn.writeline(s);
            
            %% Upload analog data
            self.conn.writeline(self.uploadAWord);
            s = sprintf(['%.6f',repmat(',%.6f',1,24),'%%'],a');
            pause(0.1);
            self.conn.writeline(s);
            
            %% Upload digital data
            self.conn.writeline(self.uploadDWord);
            s = sprintf('%d,%%',d);
            s = s(1:end-2);
            pause(0.1);
            self.conn.writeline(s);

        end
        
        function uploadDDSData(self,dds)
            %UPLOADDDSDATA Uploads compiled DDS data to a connected MOGLabs
            %ARF device.  The device to upload to must be set as the
            %RemoteControl.mog property, and you must connect to the
            %device first.
            %
            %   SELF.UPLOADDDSDATA(DDS) Uploads data in the DDS array to
            %   the MOGLabs ARF box using text-based commands only.  DDS
            %   must be a 2 element vector of struct with fields t, freq,
            %   pow, and phase.
            if isempty(self.mog)
                return
            end
            
            if isempty(self.mog.cx)
                error('Connect to MOGLabs ARF box first!');
            end
            % Create mogtable objects
            tb = mogtable(self.mog,1);
            tb(2) = mogtable(self.mog,2);
            
            % Put data into mogtable objects
            for nn = 1:numel(tb)
                tb(nn).t = dds(nn).t;
                tb(nn).freq = dds(nn).freq;
                tb(nn).pow = dds(nn).pow;
                tb(nn).phase = dds(nn).phase;
            end
            
            % Reduce instruction sizes and make sure both tables have
            % instructions at the same time
% % % % %  Use this for Gaussian Pulses
%             tb(1).reduce;
%             if sum(tb(1).sync) == 1
%                 tb(2).reduce;
%                 tb(1).reduce(tb(2).sync);
%             else
%                 tb(2).reduce(tb(1).sync);
%             end
% % % % % Use this for Square Pulses
            tb(1).reduce;
            tb(2).reduce;
% % % % %             
            % Send commands to device
            for nn = 1:numel(tb)
                self.mog.cmd('mode,%d,%s',tb(nn).channel,tb(nn).MODE);
                self.mog.cmd('table,stop,%d',tb(nn).channel);
            end
            self.mog.cmd('table,sync,1');
            numInstr = tb.upload;
            estUploadTime = numInstr*11/1000;
            if estUploadTime > (7/8*self.sq.ddsTrigDelay)
                pause(estUploadTime - self.sq.ddsTrigDelay + 1);
            end
        end
        
        function uploadBinaryDDSData(self,dds)
            %UPLOADBINARYDDSDATA Uploads compiled DDS data to a connected
            %MOGLabs ARF device using binary tables.  The device to upload
            %to must be set as the RemoteControl.mog property, and you must
            %connect to the device first.
            %
            %   SELF.UPLOADBINARYDDSDATA(DDS) Uploads data in the DDS array
            %   to the MOGLabs ARF box using binary tables.  DDS
            %   must be a 2 element vector of struct with fields t, freq,
            %   pow, and phase.
            if isempty(self.mog)
                return
            end
            
            if isempty(self.mog.cx)
                error('Connect to MOGLabs ARF box first!');
            end
            % Create mogtable objects
            tb = mogtable(self.mog,1);
            tb.pow_units = 'hex';
            tb(2) = mogtable(self.mog,2);
            tb(2).pow_units = 'hex';
            
            % Put data into mogtable objects
            for nn = 1:numel(tb)
                tb(nn).t = dds(nn).t;
                tb(nn).freq = dds(nn).freq;
                tb(nn).pow = dds(nn).pow;
                tb(nn).phase = dds(nn).phase;
            end
            
            % Send commands to device
            for nn = 1:numel(tb)
                self.mog.cmd('mode,%d,%s',tb(nn).channel,tb(nn).MODE);
                self.mog.cmd('table,stop,%d',tb(nn).channel);
            end
            tb(1).upload_binary_table;
            tb(2).upload_binary_table;
            %Re-arm the tables because stopping table 2 also stops table 1,
            %and stop command is issued in the upload process.
            for nn = 1:numel(tb)
                self.mog.cmd('table,arm,%d',tb(nn).channel);
                self.mog.cmd('table,rearm,%d,on',tb(nn).channel);
            end
            num_tries = 10;
            current_try = 1;
            while 1
                try
                    self.mog.cmd('table,sync,1');
                    break;
                catch current_exception
                    if current_try < num_tries
                        current_try = current_try + 1;
                    else
                        rethrow(current_exception);
                    end
                end
            end
            
        end
        
        function run(self,cb)
            %RUN Starts a single client run by sending the start word
            %
            %   SELF.RUN() Starts a run.
            %
            %   SELF.RUN(CB) Starts a run and sets the BytesAvailableFcn
            %   callback to the function handle CB()
            %
            self.open;
            self.conn.flush;
            if nargin > 1
                self.conn.configureCallback('terminator',@(~,~) cb());
            end
            self.conn.writeline(self.startWord);
        end %end run

        %Multicallback function 
        function multi_callback(self, callbacks)
            for k = 1:length(callbacks)
                feval(callbacks{k});
            end
        end

        function urun(self, varargin)
            %URUN Uploads current sequence and starts a run
            self.upload;
            self.run(@() self.multi_callback(varargin));
        end



        
%         function urun(self,varargin)
%             %URUN Uploads current sequence and starts a run
%             self.upload;
%             self.run(varargin{:});
%         end
        


        function loop(self,cb)
            %LOOP Starts a perpetual loop using a supplied callback
            %function
            %
            %   SELF.LOOP(CB) Starts the loop using callback function CB
            function internal_callback(~,~)
                self.read;
                cb();
                self.run;
            end
            self.conn.configureCallback('terminator',@(src,event) internal_callback);
            self.run;
        end
        
        function start(self)
            %START Starts a full run through the sequence of numRuns
            self.status = self.RUNNING;
            self.init;
            self.set;
            self.run;
        end
        
        function resume(self)
            %RESUME sets and runs a sequence
            self.status = self.RUNNING;
            self.set;
            self.run;
        end
        
        function resp(self,~,~)
            %RESP responds to the arrival a new word over TCPIP
            %   Controls the next run of the sequence, either ending it or
            %   analyzing the results and stepping forward
            s = self.read;
            if ~self.connected && strcmpi(s,self.readyWord)
                fprintf(1,'Interface connected!\n');
                self.connected = true;
                self.status = self.STOPPED;
            elseif strcmpi(s,self.readyWord) && strcmpi(self.status,self.RUNNING)
                if ~isempty(self.err) && self.err.fail()
                    % If failure, stop run
                    self.stop;
                    error('User defined error has occurred. Run stopping');
                elseif self.c.done()
                    % Analyze
                    self.analyze;
                    % Stop
                    pause(0.1);
                    self.conn.flush();
                    self.status = self.STOPPED;
                    fprintf('Run finished\n');
                    %Send a trigger to the app managing (only works with
                    %the app
                
                    appHandle = findall(0, 'Type', 'Figure', 'Name', 'Gravimeter Control Interface');
                    if isempty(appHandle)
%if no app is open nothing happens!
                    else
                        % An instance of the app already exists, so use
                        % that...
                        app = get(appHandle, 'UserData');
                    end
                    % Now we can send the trigger
                    app.triggerFromOutside();

                else
                    % Analyze
                    self.analyze;
                    % Run again
                    self.c.increment();
                    self.set;
                    self.run;
                end
            end
        end
        
        function self = init(self)
            %SET Sets the mode to INIT and calls the callback function if
            %currentRun is 1
            self.setFunc;
            if self.c.current() == 1
                self.mode = self.INIT;
                self.callback(self);
            end
        end
        
        function self = set(self)
            %SET Sets the mode to SET and calls the callback function
            self.mode = self.SET;
            self.callback(self);
        end
        
        function self = analyze(self)
            %ANALYZE Sets the mode to ANALYZE and calls the callback
            %function
            self.mode = self.ANALYZE;
            self.callback(self);
        end
        
        function r = isInit(self)
            %ISSET Returns true if the mode is INIT
            r = strcmpi(self.mode,self.INIT);
        end
        
        function r = isSet(self)
            %ISSET Returns true if the mode is SET
            r = strcmpi(self.mode,self.SET);
        end
        
        function r = isAnalyze(self)
            %ISANALYZE Returns true if the mode is ANALYZE
            r = strcmpi(self.mode,self.ANALYZE);
        end
        
        function reset(self)
            %RESET Resets currentRun to 1, data to [], mode to INIT
            self.c.reset;
            self.data = [];
            self.mode = self.INIT;
        end
        
        function r = isRunning(self)
            %ISRUNNING Returns true if the status is RUNNING
            r = strcmpi(self.status,self.RUNNING);
        end
    end %end methods

end %end classdef