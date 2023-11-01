%% Inputs
Channel = 1;
Duration = 1; % s
CentreFreq = 110; % MHz
CentrePower = 30; % dBm
CentrePhase = 0; 
RampOrPulse = 'ramp';

%% Initialise
%connect to DDS
mog = mogdevice;
mog.connect('192.168.1.100');


%clear tables
mog.cmd('MODE,1,TSB');
mog.cmd('TABLE,CLEAR,1');
mog.cmd('MODE,1,TSB');
mog.cmd('TABLE,CLEAR,2');

%% Make list of instructions
% make a pulse
if strcmpi(RampOrPulse,'ramp')
    N = 1;
    f = CentreFreq;
    P = CentrePower;
    ph = CentrePhase;
    dt = Duration*1e3;
end
% make a ramp
if strcmpi(RampOrPulse,'ramp')
    N = 10;
    f = linspace(CentreFreq - 5,CentreFreq + 5,N);
    P = ones(1,N)*CentrePower;
    ph = ones(1,N)*CentrePhase;
    dt = ones(1,N)*Duration/N*1e3;
end

%% Upload
%Display text so you know when upload starts
disp('Uploading table...')


if Channel == 1
    %Create Table to be uploaded
    for ii=1:length(N)
        % use printf notation when sending commands
        %     Channel -> Freq -> Power -> Phase -> Duration
        mog.cmd('TABLE,APPEND,1,%f MHz,%f dBm,%f,%f ms',f(ii),P(ii),ph(ii),dt(ii));
    end
    mog.cmd('TABLE,APPEND,1,110 MHz,0 dBm,0,1us');

    % arm dds
    mog.cmd('TABLE,ARM,1');    
else
    for ii=1:length(N)
        % use printf notation when sending commands
        %     Channel -> Freq -> Power -> Phase -> Duration
        mog.cmd('TABLE,APPEND,1,%f MHz,%f dBm,%f,%f ms',f(ii),P(ii),ph(ii),dt(ii));
    end
    mog.cmd('TABLE,APPEND,1,110 MHz,0 dBm,0,1us');

    % arm dds
    mog.cmd('TABLE,ARM,2');       
end
    
%Display text so you know when the upload/arming is done
disp('Done')



%Run through the table
    % run channel 1 commands
mog.cmd('TABLE,START,1') 

%Disconnect
% delete(mog);
clear Channel Duration CentreFreq CentrePhase CentrePower RampOrPulse
clear N f P ph dt

%% old
% 
% %connect to DDS
% mog = mogdevice;
% mog.connect('192.168.1.100');
% 
% %Ramp details
% N = 10;
% fRamp = linspace(105,115,N);
% 
% %clear table
% mog.cmd('MODE,1,TSB');
% mog.cmd('TABLE,CLEAR,1');
% 
% %Display text so you know when upload starts
% disp('Uploading table...')
% 
% %Create Table to be uploaded
% for i=1:length(fRamp)
%     % we can use printf notation when sending commands
%     mog.cmd('TABLE,APPEND,1,%f MHz,30 dBm,0,500 ms',fRamp(i));
% end
%     mog.cmd('TABLE,APPEND,1,115 MHz,0 dBm,0,1us');
%     
% %Display text so you know when the upload is done
% disp('Done')
% 
% % arm dds
% mog.cmd('TABLE,ARM,1');
% %Run through the table
% % mog.cmd('TABLE,START,1')
% 
% %Disconnect
% % delete(mog);
% 
