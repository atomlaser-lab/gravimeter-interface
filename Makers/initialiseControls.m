% Function for controls initialisation of the Gravy
function initialiseControls()

% Create remote control instance
r = RemoteControl;

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
disp(['Moglabs DDS connected: ', mogipAddress])

%% MKU device initialization
% Define the mku MAC address
mkumacAddress = 'e8-dB-84-dA-b5-02';

% Get the IP address associated with the mku MAC address
mkuipAddress = MACtoIPFinder(mkumacAddress,'no ping');

% Store in IPaddresses struct
IPaddresses.MKU = mkuipAddress;

% Create mku instance and connect it
m = mku(mkuipAddress);
m.writeList((const.f_Rb_groundHFS/1e6 + [-315e-3+4e-3,0.069e-3])/2*1e6); %values subject to change day to day and dependng on mag field
% Display a success message for mku connection
disp(['Microwave Synthetizer mku connected. IP address: ', mkuipAddress])

%% Opt Initialisation
opt = SequenceOptions;
r.devices.opt = opt;

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
assignin('base', 'r', r);
assignin('base', 'mog', mog);
assignin('base', 'm', m);
assignin('base', 'opt', opt);
assignin('base', 'IPaddresses', IPaddresses);

disp('IP addresses have been saved in the workspace');

%% Abs_Analysis_GUI initialization
% If 'Abs Analysis Parameters Loader' window exists, close it
existingWindow1 = findobj('Type', 'figure', 'Name', 'Abs Analysis Parameters Loader', 'Visible', 'on');
if ~isempty(existingWindow1)
    close(existingWindow1);
end

% If 'Abs Analysis GUI Figure' window exists, close it
existingWindow2 = findobj('Type', 'figure', 'Name', 'Abs Analysis GUI Figure', 'Visible', 'on');
if ~isempty(existingWindow2)
    close(existingWindow2);
end

% If AbsrotionImage class object exists in workspace, delete it
if evalin('base', 'exist(''AbsrotionImage'', ''var'')')
    evalin('base', 'clear AbsrotionImage');
end

% If Abs_Analysis_parameters struct exists in workspace, delete it
if evalin('base', 'exist(''Abs_Analysis_parameters'', ''var'')')
    evalin('base', 'clear Abs_Analysis_parameters');
end

% Open/reopen 
% Abs_Analysis_GUI;

% disp('Abs Anaclysis for Absorption Imaging Initialized')

% Display the remote control and options
disp(r)
disp(opt)
end
