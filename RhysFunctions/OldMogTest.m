% old

%connect to DDS
mog = mogdevice;
mog.connect('192.168.1.100');

%Ramp details
N = 10;
fRamp = linspace(105,115,N);

%clear table
mog.cmd('MODE,1,TSB');
mog.cmd('TABLE,CLEAR,1');

%Display text so you know when upload starts
disp('Uploading table...')

%Create Table to be uploaded
for i=1:length(fRamp)
    % we can use printf notation when sending commands
    mog.cmd('TABLE,APPEND,1,%f MHz,30 dBm,0,500 ms',fRamp(i));
end
%     mog.cmd('TABLE,APPEND,1,115 MHz,0 dBm,0,1us');
    
%Display text so you know when the upload is done
disp('Done')

% arm dds
mog.cmd('TABLE,ARM,1');
%Run through the table
% mog.cmd('TABLE,START,1')

%Disconnect
% delete(mog);

