% Function for controls initialisation of the Gravy
function initialiseControls(DDS)
closeNoPrompt(matlab.desktop.editor.getAll);
% open makeBEC_RamanInterferometer.m
% open Callback_Rhys_2State
% open Abs_Analysis_DualState_RT.m
% open Callback_Rhys_2State_MW.m

open makeBEC_Bragg_NewCloud.m
open Callback_Rhys.m
open Abs_Analysis_GUI.m
open Callback_Rhys_2ROI_MW.m

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
if strcmpi(DDS,'Bragg')
    % mog.connect(mogipAddress);
    mog.connect('192.168.1.4');
    r.mog = mog;
    % Display a success message for moglabs connection
    % disp(['Moglabs DDS connected: ', mogipAddress])
    disp(['Moglabs DDS connected: ', '192.168.1.4'])
else
    mog.connect('192.168.1.4');
    r.mog = mog;
end
%% MKU device initialization

% % % Just enter the address to reduce start up time
m = mku('192.168.1.10');
m.writeList((const.f_Rb_groundHFS - 315e3 - 6.9e3)/2,(const.f_Rb_groundHFS)/2);
r.devices.mku = m;
disp(['Microwave Synthetizer mku connected. IP address: ', '192.168.1.10'])

%% Opt Initialisation
opt = SequenceOptions;
r.devices.opt = opt;

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
if strcmpi(DDS,'Bragg')
    Abs_Analysis_GUI;
    set(figure(2),'WindowStyle','Docked');
else
    close all
%     Abs_Analysis_DualState_RT('last',1);
    Abs_Analysis_GUI('last',1);    
    set(figure(1),'WindowStyle','Docked');
    set(figure(2),'WindowStyle','Docked');
end
% Display the remote control and options
disp(r)
disp(opt)

%% set callback
r.callback = @Callback_Rhys_2ROI_MW;
r.makerCallback = @makeBEC_Bragg_NewCloud;
opt.mw.enable(1) = 1;
opt.mw.enable_sg = 1;
opt.mw.analyze(1) = 1;
end
