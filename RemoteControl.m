classdef RemoteControl < handle
    properties
        %% TCPIP properties
        conn        %TCPIP connection
        connected   %Is LabVIEW client connected?
        %% Sequence properties
        status      %Current status of run: RUNNING or STOPPED
        currentRun  %Current run number
        numRuns     %Total number of runs
        %% Data properties
        mode        %Mode of callback function: SET, ANALYZE, or INIT
        devices     %Structure listing MATLAB devices used in callback
        data        %Data structure to use in callback function
        callback    %Callback function, takes argument of Rebeka selfect
    end
    
    properties(SetAccess = protected)
        remoteAddress = 'localhost';    %Connect to local host
        remotePort = 6666;            %Remote port to use
        waitTime = 0.1;               %Wait time for querying BytesAvailable and between successive writes

    end %end constant properties
    
    properties(Constant, Hidden=true)
        readyWord = 'ready';          %Word indicating that client is ready
        startWord = 'start';          %Word telling host to start
        endWord = 'end';              %Word telling host to stop TCP loop
        uploadDWord = 'uploadD';      %Word telling host to upload digital (uint32) data
        uploadAWord = 'uploadA';      %Word telling host to upload analog (float) data

        SET = 'set/check';
        ANALYZE = 'analyze';
        INIT = 'init';
        
        RUNNING = 'running';
        STOPPED = 'stopped';
    end        
    
    methods
        function self = RemoteControl
            self.connected = false;
            self.mode = self.INIT;
            self.status = self.STOPPED;
            self.reset;
        end %end constructor
        
        function open(self)
            %OPEN Opens a tcpip port           
            %open Creates and opens a TCP conn.  Waits for ready word
            if isempty(self.conn) || ~isvalid(self.conn)
                self.conn = tcpip(self.remoteAddress,self.remotePort,'networkrole','client');
                self.conn.Terminator = 'CR/LF';
                self.conn.BytesAvailableFcn = @(src,event) self.resp(src,event);
                self.connected = false;
                self.conn.OutputBufferSize = 2^24;
            end
            
            if strcmpi(self.conn.Status,'closed')
                fprintf(1,'Attempting connection...\n');
                fopen(self.conn);
                fprintf(1,'Connection successful!\n');
                self.connected = true;
            end
        end %end open
        
        function setFunc(self)
            self.conn.BytesAvailableFcn = @(src,event) self.resp(src,event);
        end
        
        function r = findTCPPort(self)
            %FINDTCPPORT Finds all existing TCP ports
            instr = instrfindall;
            for nn = 1:numel(instr)
                if isa(instr(nn),'tcpip') && strcmpi(instr(nn).RemotePort,self.remotePort) && strcmpi(instr(nn).RemoteHost,self.remoteAddress)
                    r = instr(nn);
                    return;
                end
            end
            r = false;       
        end
        
        function r = read(self)
            %READ Reads available data from TCP connection
            r = fgetl(self.conn);
        end %end read
        
        function r = waitForReady(self)
            %WAITFOREADY Returns true when the readyWord appears
            while self.conn.BytesAvailable <= length(self.readyWord)
                pause(self.waitTime); 
            end
            r = strcmpi(self.read,self.readyWord);
        end %end waitForReady
        
        function stop(self)
            %STOP Releases client from remote control and closes TCP conn
            if strcmpi(self.conn.Status,'open')
%                 fprintf(self.conn,'%s',self.endWord);
                fclose(self.conn);
            end
            delete(self.conn);
            fprintf(1,'Remote control session terminated\n');
            self.connected = false;
            self.status = self.STOPPED;
        end %end fclose
        
        function upload(self,data)
            %UPLOAD uploads data to host
            %
            %   upload(r,data) with r the RemoteControl object and data a
            %   2D array with times in the first column, a 32 bit digital
            %   value in the second column, and 24 analog values in the rest
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
            
            self.open;
            
            fprintf(self.conn,'%s\n',self.uploadAWord);
            s = sprintf(['%.6f',repmat(',%.6f',1,24),'%%'],a');
            pause(0.1);
            fprintf(self.conn,s);
            
            fprintf(self.conn,'%s\n',self.uploadDWord);
            s = sprintf('%d,%%',d);
            s = s(1:end-3);
            pause(0.1);
            fprintf(self.conn,s);
            
        end
        
        function run(self)
            %RUN Starts a single client run by sending the start word
            self.open;
            fprintf(self.conn,'%s\n',self.startWord);
        end %end run
        
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
                if self.currentRun == self.numRuns
                    % Analyze
                    self.analyze;
                    % Stop
                    self.stop;
                else
                    % Analyze
                    self.analyze;
                    % Run again
                    self.currentRun = self.currentRun + 1;
                    self.set;
                    self.run;
                end
            end
        end
        
        function self = init(self)
            %SET Sets the mode to INIT and calls the callback function if
            %currentRun is 1
            if self.currentRun==1
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
            self.currentRun = 1;
            self.data = [];
            self.mode = self.INIT;
        end
        
    end %end methods

end %end classdef