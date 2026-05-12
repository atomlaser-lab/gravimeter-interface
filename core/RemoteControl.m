classdef RemoteControl < handle
    %REMOTECONTROL Class for handling remote connections to the Quantum
    %Sensors group's LabVIEW control interface.
    properties
        % TCPIP properties
        conn            %TCPIP connection
        connected       %Is LabVIEW client connected?
        % Sequence properties
        status          %Current status of run: RUNNING or STOPPED
        sq              %Sequence object representing current sequence
        makerCallback   %Callback function for creating a TimingSequence object
        % Data properties
        mode            %Mode of callback function: SET, ANALYZE, or INIT
        devices         %Structure listing MATLAB devices used in callback
        data            %Data structure to use in callback function
        callback        %Callback function, takes argument of Rebeka object
    end
    
    properties(SetAccess = protected)
        remoteAddress = 'localhost';        %Connect to local host
        remotePort = 6666;                  %Remote port to use
    end

    properties(Access = protected)
        run_callback                        %Callback to use for run() and loop()
        wait_for_image                      %Flag to indicate that we need to wait for imaging to be done
        maker_copy                          %Copy of the maker callback function with arguments
    end
    
    properties(SetAccess = immutable)
        c                                   %Rollover counter object, keeps track of runs
    end
    
    properties(Constant, Hidden=true)
        CMD_READY = 'ready';                %Word indicating that client is ready
        CMD_START = 'start';                %Word telling host to start
        CMD_END = 'end';                    %Word telling host to stop TCP loop
        CMD_UPLOAD_DIGITAL = 'uploadD';     %Word telling host to upload digital (uint32) data
        CMD_UPLOAD_ANALOG = 'uploadA';      %Word telling host to upload analog (float) data
        CMD_CAM_DELAY = 'camDelay';         %Word telling host to store camera acquisition delay
        CMD_WAIT_FOR_IMAGE = 'waitForImage';%Word telling host to wait for image before sending end word

        CAM_STATUS_FMT = 'status: %d, image: %d';   %Format for camera status information
        CAM_STATUS_NO_ERR = 0;
        CAM_STATUS_TIMEOUT = -1;

        MODE_SET = 'set/check';             %Indicates that the callback mode is to set parameters
        MODE_ANALYZE = 'analyze';           %Indicates that the callback mode is to analyze data
        MODE_INIT = 'init';                 %Indicates that the callback mode is to initialize the run
        
        STATUS_AUTO = 'auto';               %Indicates that an automated sequence is running
        STATUS_LOOP = 'loop';               %Indicates that an automated loop is running
        STATUS_STOPPED = 'stopped';         %Indicates that an automated sequence is not running

        MAKER_STORAGE_DIRECTORY = 'D:\run-files';           %Directory for storing copies of maker functions
        MAKER_FILENAME_FORMAT = '%s_Image%d.m';             %String format for maker function copies
        OPTIONS_FILENAME_FORMAT = 'options_Image%d.mat';    %String format for options filenames
    end    
    
    events
        sequenceChanged                     %Event for notifying that a sequence has changed
    end
    
    methods
        function self = RemoteControl(varargin)
            %REMOTECONTROL Constructs a RemoteControl object
            %
            %   SELF = REMOTECONTROL(ADDRESS) Constructs an object that
            %   will connect to the control VI at ADDRESS on the default
            %   port
            %
            %   SELF = REMOTECONTROL(ADDRESS,PORT) object will connect
            %   using given address and port.
            self.setRemoteProperties(varargin{:});
            self.connected = false;
            self.mode = self.MODE_INIT;
            self.status = self.STATUS_STOPPED;
            self.makerCallback = @makeSequenceGinger;
            self.c = RolloverCounter();
            self.reset;
        end
        
        function self = setRemoteProperties(self,varargin)
            %SETREMOTEPROPERTIES Sets the remote address and port
            %
            %   SELF = SETREMOTEPROPERTIES(ADDRESS,PORT) sets the remote
            %   address and (optional) port
            if numel(varargin) >= 1
                self.remoteAddress = varargin{1};
            end
            if numel(varargin) >= 2
                self.remotePort = varargin{2};
            end
        end
        
        function open(self)
            %OPEN Opens a tcpip port
            %
            %   OPEN() Creates and opens a TCP connection, sets connected
            %   property to TRUE
            if isempty(self.conn)
                fprintf(1,'Attempting connection...\n');
                self.conn = tcpclient(self.remoteAddress,self.remotePort);
                self.conn.configureTerminator('CR/LF');
                self.conn.configureCallback('terminator',@(src,event) self.resp(src,event));
                self.conn.ErrorOccurredFcn = @(src,event) RemoteControl.error_handler(src,event);
                R = version('-release');
                release_year = regexp(R,'\d+','match');
                release_year = str2double(release_year{1});
                if contains(R,'2022b') || (release_year > 2022)
                    self.conn.InputBufferSize = 2^20;
                    self.conn.OutputBufferSize = 2^20;
                end
                fprintf(1,'Connection successful!\n');
                self.connected = true;
            end
        end %end open
        
        function setFunc(self)
            %SETFUNC Sets the BytesAvailableFcn to self.resp()
            self.open;
            self.conn.configureCallback('terminator',@(src,event) self.resp(src,event))
        end

        function loop(self,cb)
            %LOOP Starts a perpetual loop
            %
            %   SELF = SELF.LOOP(CB) Creates a perpetual loop that calls
            %   the function handle CB on every run.
            self.open;
            if nargin < 2
                self.run_callback = cb;
            end
            self.status = self.STATUS_LOOP;
            self.setFunc;
            self.run;
        end
        
        function r = read(self)
            %READ Reads available data from TCP connection
            r = self.conn.readline();
        end
        
        function stop(self)
            %STOP Releases client from remote control and closes TCP
            %connection
            self.conn = [];
            fprintf(1,'Remote control session terminated\n');
            self.connected = false;
            self.status = self.STATUS_STOPPED;
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
            %   SELF = MAKE(SELF,VARARGIN) runs
            %   SELF.MAKERCALLBACK(VARARGIN{:}) and stores the resulting
            %   sequence in SELF.SQ.
            %
            %   Notifies listeners that the "sequenceChanged" event has
            %   occurred.

            if isempty(self.makerCallback) || ~isa(self.makerCallback,'function_handle')
                error('Provide a valid sequence creation function to makerCallback!');
            end
            self.sq = self.makerCallback(varargin{:});
            notify(self,'sequenceChanged');
            [self.maker_copy.maker,self.maker_copy.opt] = copy_sequence(self.makerCallback,varargin{:});
        end
        
        function self = upload(self,data)
            %UPLOAD uploads data to host
            %
            %   SELF = UPLOAD uploads data to control interface using the 
            %   current sequence stored in the SELF.SQ field
            %
            %   SELF = UPLOAD(DATA) with uploads data structure DATA.  DATA
            %   must be a 2D array with times in the first column, a 32 bit
            %   digital value in the second column, and 24 analog values in
            %   the rest
            if nargin < 2
                data = self.sq.compile;
                self.wait_for_image = data.waitForImage;
            else
                self.wait_for_image = data.waitForImage;
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
            
            %% Check TA status, turn on TA if necessary
            ldd = self.get_devices('mog');
            for nn = 1:numel(ldd)
                ldd{nn}.turn_on;
            end
            %% Upload DDS data
            self.uploadDDSData(data.dds);
            
            %% Open connection with LabVIEW VI and set options
            self.open;
            % This does the camera delay value
            self.conn.writeline(self.CMD_CAM_DELAY);
            s = sprintf('%.1f',data.camDelay);
            pause(0.1);
            self.conn.writeline(s);
            % This does the waitForImage flag
            self.conn.writeline(self.CMD_WAIT_FOR_IMAGE);
            s = sprintf('%.0f',data.waitForImage);
            pause(0.1);
            self.conn.writeline(s);
            
            %% Upload analog data
            self.conn.writeline(self.CMD_UPLOAD_ANALOG);
            s = sprintf(['%.6f',repmat(',%.6f',1,24),'%%'],a');
            pause(0.1);
            self.conn.writeline(s);
            
            %% Upload digital data
            self.conn.writeline(self.CMD_UPLOAD_DIGITAL);
            s = sprintf('%d,%%',d);
            s = s(1:end-2);
            pause(0.1);
            self.conn.writeline(s);

        end
        
        function uploadDDSData(self,dds)
            %UPLOADDDSDATA Uploads the DDS data via the MOGLABS interface
            %
            %   UPLOADDDSDATA(DDS) uploads DDS data stored in DDS
            mog = self.get_devices('mogrf');
            if isempty(mog)
                return
            elseif numel(mog) > 1
                error('More than one MOGRF object is not supported when uploading DDS data!');
            else
                mog = mog{1};
            end
            
            if isempty(mog.cx)
                error('Connect to MOGLabs ARF box first!');
            end
            % Create mogtable objects
            tb = mogtable(mog,1);
            tb(2) = mogtable(mog,2);
            tb(1).pow_units = 'hex';
            tb(2).pow_units = 'hex';
            
            % Put data into mogtable objects
            for nn = 1:numel(tb)
                tb(nn).t = dds(nn).t;
                tb(nn).freq = dds(nn).freq;
                tb(nn).pow = dds(nn).pow;
                tb(nn).phase = dds(nn).phase;
            end
            
            % Reduce instruction sizes and make sure both tables have
            % instructions at the same time
            if ~strcmpi(tb(1).pow_units,'hex') && strcmpi(tb(2).pow_units,'hex') 
                tb(1).reduce;
                if sum(tb(1).sync) == 1
                    tb(2).reduce;
                    tb(1).reduce(tb(2).sync);
                else
                    tb(2).reduce(tb(1).sync);
                end
            end
            
            % Send commands to device
            for nn = 1:numel(tb)
                mog.cmd('mode,%d,%s',tb(nn).channel,tb(nn).MODE);
                mog.cmd('table,stop,%d',tb(nn).channel);
            end
            num_tries = 10;
            current_try = 1;
            while 1
                try
                    mog.cmd('table,sync,1');
                    break;
                catch err
                    if current_try < num_tries
                        current_try = current_try + 1;
                    else
                        rethrow(err);
                    end
                end
            end
            numInstr = tb.upload;
            estUploadTime = numInstr*11/3280;
            if estUploadTime > (7/8*self.sq.ddsTrigDelay)
                pause(estUploadTime - self.sq.ddsTrigDelay + 1);
            end
        end
        
        function run(self,cb)
            %RUN Starts a single client run by sending the start word
            self.open;
            self.conn.flush;
            if nargin > 1
                self.run_callback = cb;
            end
            self.conn.writeline(self.CMD_START);
            %====================================SAM HACK HERE===================================
%             addpath('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\SDK');
%             pause(18)
%             Blink_SDK_Sam_v1;
            %==================================== end ==========================================
        end %end run

        function urun(self,varargin)
            %URUN Uploads current sequence and starts a run
            self.upload;
            self.run(varargin{:});
        end
        
        function start(self)
            %START Starts a full run through the sequence of numRuns
            self.status = self.STATUS_AUTO;
            self.init;
            self.set;
            self.run;
        end
        
        function resume(self)
            %RESUME sets and runs a sequence
            self.status = self.STATUS_AUTO;
            self.set;
            self.run;
        end
        
        function resp(self,~,~)
            %RESP responds to the arrival a new word over TCPIP
            %   Controls the next run of the sequence, either ending it or
            %   analyzing the results and stepping forward
            
            %
            % First, we grab the camera error information and image number,
            % if present
            %

            s = self.read;
            if self.wait_for_image && ~strcmpi(s,self.CMD_READY)
                % This executes if we need to wait for image acquisition to
                % complete, and what is sent by the control VI is not the
                % ready word
                r = sscanf(s,self.CAM_STATUS_FMT);
                % Check camera error codes
                if r(1) == self.CAM_STATUS_TIMEOUT
                    error('Camera acquisition timed out');
                elseif r(1) ~= self.CAM_STATUS_NO_ERR
                    error('Unknown camera error');
                end
                % If no error, grab the image number
                image_number = r(2);
                % Save a copy of the maker function with the associated
                % image number
                fid = fopen(self.get_maker_filename(func2str(self.makerCallback),image_number),'w');
                fprintf(fid,'%s',self.maker_copy.maker);
                fclose(fid);
                % Save a copy of the standard sequence options as a MAT
                % file
                opt = self.maker_copy.opt;
                save(self.get_options_filename(image_number),'opt');
            elseif strcmpi(s,self.CMD_READY)
                % If the Control VI sends CMD_READY, execute the
                % appropriate callback function
                if strcmpi(self.status,self.STATUS_AUTO)
                    % Status is AUTO only for automated sequences
                    if self.c.done()
                        % Analyze
                        self.analyze;
                        % Stop
                        pause(0.1);
                        self.conn.flush;
                        self.status = self.STATUS_STOPPED;
                        fprintf(1,'Run finished\n');
                    else
                        % Analyze
                        self.analyze;
                        % Run again
                        self.c.increment();
                        self.set;
                        self.run;
                    end
                elseif ~isempty(self.run_callback) && isa(self.run_callback,'function_handle')
                    % If not an automated sequence, and a run_callback
                    % function is present, run that callback
                    self.run_callback();
                end
                % If we are looping, then start another run. This occurs
                % even if no run_callback is present
                if strcmpi(self.status,self.STATUS_LOOP)
                    self.run;
                end
            end
            % Call resp() again if there are still data to be read
            if self.conn.BytesAvailable
                self.resp();
            end
        end
        
        function self = init(self)
            %SET Sets the mode to INIT and calls the callback function if
            %currentRun is 1
            self.setFunc;
            if self.c.current() == 1
                self.mode = self.MODE_INIT;
                self.callback(self);
            end
        end
        
        function self = set(self)
            %SET Sets the mode to SET and calls the callback function
            self.mode = self.MODE_SET;
            self.callback(self);
        end
        
        function self = analyze(self)
            %ANALYZE Sets the mode to ANALYZE and calls the callback
            %function
            self.mode = self.MODE_ANALYZE;
            self.callback(self);
        end
        
        function r = isInit(self)
            %ISSET Returns true if the mode is INIT
            r = strcmpi(self.mode,self.MODE_INIT);
        end
        
        function r = isSet(self)
            %ISSET Returns true if the mode is SET
            r = strcmpi(self.mode,self.MODE_SET);
        end
        
        function r = isAnalyze(self)
            %ISANALYZE Returns true if the mode is ANALYZE
            r = strcmpi(self.mode,self.MODE_ANALYZE);
        end
        
        function reset(self)
            %RESET Resets currentRun to 1, data to [], mode to INIT
            self.c.reset;
            self.data = [];
            self.mode = self.MODE_INIT;
        end

        function particular_devices = get_devices(self,device_class)
            %GET_DEVICE Returns all instances of a particular device class
            all_names = fieldnames(self.devices);
            particular_devices = {};
            for nn = 1:numel(all_names)
                current_device = self.devices.(all_names{nn});
                if isa(current_device,device_class)
                    particular_devices{end + 1} = current_device; %#ok<AGROW> 
                end
            end
        end
        
    end

    methods(Static)
        function fname = get_maker_filename(maker_filename,image_number)
            fname = sprintf(RemoteControl.MAKER_FILENAME_FORMAT,maker_filename,image_number);
            fname = fullfile(RemoteControl.MAKER_STORAGE_DIRECTORY,fname);
        end

        function fname = get_options_filename(image_number)
            fname = fullfile(RemoteControl.MAKER_STORAGE_DIRECTORY,sprintf(RemoteControl.OPTIONS_FILENAME_FORMAT,image_number));
        end

        function error_handler(src,event)
            src
            event
        end
    end

end