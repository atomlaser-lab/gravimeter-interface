% Function for debugging ports when IP address changes abruptly
function DebugPorts()
clear mog
clear m
clear IPaddresses
r.devices = [];
r.mog = [];
% Create IPaddresses struct
IPaddresses = struct();

%% Moglabs device initialization
% Define the moglabs MAC address
mogmacAddress = '70-b3-d5-84-aa-3a';

% Get the IP address associated with the moglabs MAC address
mogipAddress = MACtoIPFinder(mogmacAddress);

% Store in IPaddresses struct
IPaddresses.MogDDS = mogipAddress;

% Create moglabs device instance and connect it
mog = mogdevice;
mog.connect(mogipAddress);
r.mog = mog;
% Display a success message for moglabs connection
disp(['Moglabs DDS re-connected: ', mogipAddress])

%% MKU device initialization
% Define the mku MAC address
mkumacAddress = 'e8-dB-84-dA-b5-02';

% Get the IP address associated with the mku MAC address
mkuipAddress = MACtoIPFinder(mkumacAddress,'no ping');

% Store in IPaddresses struct
IPaddresses.MKU = mkuipAddress;

% Create mku instance and connect it
m = mku(mkuipAddress);

% Display a success message for mku connection
disp(['Microwave Synthetizer mku re-connected. IP address: ', mkuipAddress])



%% Find other devices
% Find the IP of the MTS setup redpitaya (currently rp-f082a6)
MTSmacaddress = '00-26-32-f0-82-a6';
MTSipAddress = MACtoIPFinder(MTSmacaddress,'no ping');

% Store in IPaddresses struct
IPaddresses.MTS = MTSipAddress;

% Find the IP of the FMI setup redpitaya (currently rp-f0919a)
FMImacaddress = '00-26-32-f0-91-9a';
FMIipAddress = MACtoIPFinder(FMImacaddress,'no ping');

% Store in IPaddresses struct
IPaddresses.FMI = FMIipAddress;

%% Assign the variables r, mog, and m in the base workspace
assignin('base', 'mog', mog);
assignin('base', 'm', m);
assignin('base', 'IPaddresses', IPaddresses);

disp('IP addresses have been saved in the workspace');
end
