% Function for controls initialisation of the Gravy
function initialiseControlsYosri()
open CloseMZ_IntensityNoise.mlx
%Check the input argument
narginchk(0,1)

if   nargin<1
    DDS = 'Raman';
end

% Create remote control instance
r = RemoteControl;

% Create IPaddresses struct
IPaddresses = struct();

%% Moglabs device initialization

% Create moglabs device instance and connect it
    mog = mogdevice;
    mog.connect('192.168.1.4');
    r.mog = mog;
    % Display a success message for moglabs connection
    % disp(['Moglabs DDS connected: ', mogipAddress])
    disp(['Moglabs DDS connected: ', '192.168.1.4'])

%% MKU device initialization

% % % Just enter the address to reduce start up time
m = mku('192.168.1.10');
f1 = (const.f_Rb_groundHFS/1e6 - 315e-3 + (17-6.3)*1e-3)/2*1e6;
m.writeList(f1,(const.f_Rb_groundHFS)/2);
r.devices.mku = m;
disp(['Microwave Synthetizer mku connected. IP address: ', '192.168.1.10'])

%% Opt Initialisation
opt = SequenceOptions;
r.devices.opt = opt;
opt.MOT_LoadTime =12;
opt.detuning = 0;
opt.TwoStateImaging = 0;
opt.dipoles = 1.6;
opt.tof = 216.5e-3;
r.makerCallback = @CloseMZ_IntensityNoise;

%% Find other devices
% % Find the IP of the MTS setup redpitaya (currently rp-f082a6)
% MTSmacaddress = '00-26-32-f0-82-a6';
% MTSipAddress = MACtoIPFinder(MTSmacaddress,'no ping');
%
% % Store in IPaddresses struct
% IPaddresses.MTS = MTSipAddress;
%
% % Find the IP of the FMI setup redpitaya (currently rp-f0919a)
% FMImacaddress = '00-26-32-f0-91-9a';
% FMIipAddress = MACtoIPFinder(FMImacaddress,'no ping');
%
% % Store in IPaddresses struct
% IPaddresses.FMI = FMIipAddress;

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
    Abs_Analysis_GUI;
    set(figure(2),'WindowStyle','Docked');
% Display the remote control and options
disp(r)
disp(opt)

end
