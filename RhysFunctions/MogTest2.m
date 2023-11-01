%connect to DDS
dev = mogdevice();
dev.connect('192.168.1.100');

%Ramp details
N = 1000;
fRamp = linspace(80,300,N);

%clear table
dev.cmd('MODE,1,TSB');
dev.cmd('TABLE,CLEAR,1');

%Display text so you know when upload starts
disp('Uploading table...')

%Create Table to be uploaded
for i=1:length(fRamp)
    % we can use printf notation when sending commands
    dev.cmd('TABLE,APPEND,1,%f MHz,1 dBm,0,500 us',fRamp(i));
end
    dev.cmd('TABLE,APPEND,1,300 MHz,1 dBm,0,1us')
    
%Display text so you know when the upload is done
disp('Done')

%Run through the table
dev.cmd('TABLE,ARM,1')
dev.cmd('TABLE,START,1')

%Disconnect
delete(dev);
