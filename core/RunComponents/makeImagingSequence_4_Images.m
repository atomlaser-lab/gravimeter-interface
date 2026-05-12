function makeImagingSequence_4_Images(sq,varargin)

sq.waitForImage = true;
%
% Define default parameters
%
pulseTime = 30e-6;
pulse_delay = 15e-6;
repumpTime = 100e-6;
repumpDelay = 5e-3;
repumpShutterDelay = 5e-3;
camTime = 100e-6;
cycleTime = 40e-3;
repumpFreq = 0;
repumpAmplitude = 1;
imgFreq = 8.5;
imgAmplitude = 11;
imaging_field = 1;
image_type = 'horizontal';
%
% Parse input arguments as name/value pairs
%
if mod(numel(varargin),2) ~= 0
    error('Input arguments must be in name/value pairs');
else
    for nn = 1:2:numel(varargin)
        p = lower(varargin{nn});
        v = varargin{nn+1};
        switch p
            case 'tof'
                tof = v;
            case 'tof2'
                tof2 = v;
            case 'pulse time'
                pulseTime = v;
            case 'pulse delay'
                pulse_delay = v;
            case 'repump time'
                repumpTime = v;
            case 'repump delay'
                repumpDelay = v;
            case 'cycle time'
                cycleTime = v;
            case 'cam time'
                camTime = v;
            case 'repump freq'
                repumpFreq = v;
            case 'repump amp'
                repumpAmplitude = v;
            case 'imaging freq'
                imgFreq = v;
            case 'imaging amplitude'
                imgAmplitude = v;
            case 'repump shutter delay'
                repumpShutterDelay = v;
            case 'imaging_field'
                imaging_field = v;
            case 'image type'
                image_type = v;
            otherwise
                error('Unsupported option %s',p);
        end
    end
end

%save curennt time as drop time again
timeAtDrop = sq.time;

% Set imaging parameters BEFORE you take the image
if strcmpi(image_type , 'vertical')
    sq.find('CD0 Fast').after(tof,0.8); %zero mag field
else
sq.find('CD0 Fast').after(tof+1e-3,0); %zero mag field
sq.find('MOT bias coil').set(imaging_field); %turn on the imaging coil (to align the axis of atoms)
% sq.find('MOT bias coil').after(tof,imaging_field); %turn on the imaging coil (to align the axis of atoms)
sq.find('MOT bias').after(tof-2e-3,1); %ttl on imaging coil %tof-12.5e-3
end

%
% Preamble - set the imaging frequency
%
sq.find('87 imag freq').set(imgFreq);
sq.find('87 imag amp').set(imgAmplitude);
%
% Set camera type
%
if strcmpi(image_type,'horizontal')
    cam_trig = '87 cam trig';
    img_ch = '87 imag';
elseif strcmpi(image_type , 'vertical')
    cam_trig = 'ND cam trig'; %vertical cam trig
    img_ch = 'ND imag';
elseif strcmpi(image_type,'85')
    cam_trig = 'ND cam trig';
elseif strcmpi(image_type,'MOT')
    cam_trig = '87 cam trig';
else
    warning('incompatible cam trig input')
end

imageF2_time = tof;
imageF1_time = tof+tof2;
repump_time = tof+tof2-repumpTime-repumpDelay;
%repump_time = tof-repumpTime-repumpDelay;

%
% Imaging beam and camera trigger for image with atoms in F = 2 state
%
sq.anchor(timeAtDrop);
sq.find(img_ch).after(imageF2_time,1).after(pulseTime,0); %Turn on after TOF, then turn off after pulse time
sq.find(cam_trig).after(imageF2_time - pulse_delay,1).after(camTime,0);    %Turn on after TOF, then turn off after camera time

% %blow away F=2 atoms (set up on 28/05/2024)
% Tblow = 10e-6;
% sq.find('3DMOT').after(tof+1e-3,1);
% sq.delay(Tblow);
% sq.find('3DMOT').after(tof+1e-3+Tblow,0);

%old way to blow away F=2 atoms (up until 28/05/2024)
sq.find(img_ch).after(1e-3,1).after(1e-3,0); %get rid of f=2 atoms

% sq.waitFromLatest(cycleTime);                       %Delay
%
% Set repump values to pump F = 1 atoms into F = 2
%
sq.anchor(timeAtDrop);
sq.find('87 repump freq').after(repump_time,repumpFreq);
sq.find('87 repump amp').after(repump_time,repumpAmplitude);
%Turn on the repump TTL and the shutter
sq.find('87 repump').after(repump_time,1).after(repumpTime,0);
sq.find('Repump shutter').after(tof + pulseTime,1);
%
% Imaging beam and camera trigger for image with atoms in F=1 (moved to F=2 by repump)
%
sq.anchor(timeAtDrop);
sq.find(img_ch).after(imageF1_time,1).after(pulseTime,0); %Turn on after TOF, then turn off after pulse time
sq.find(cam_trig).after(imageF1_time - pulse_delay,1).after(camTime,0);    %Turn on after TOF, then turn off after camera time
sq.waitFromLatest(cycleTime);                       %Delay
%
% Take image without atoms
%
sq.find(img_ch).set(1).after(pulseTime,0);       %Turn on after TOF, then turn off after pulse time
sq.find(cam_trig).before(pulse_delay,1).after(camTime,0);          %Turn on after TOF, then turn off after pulse time
sq.waitFromLatest(cycleTime);                       %Delay
sq.find('Repump shutter').before(50e-6,0);    %Turn off fiber switch
%
% Take a dark image
%
sq.find(cam_trig).set(1).after(camTime,0);   %Turn on after TOF, then turn off after camera time
sq.anchor(sq.latest);   %Re-anchor the sequence to the latest value
%
% Need a last instruction so that the run ends properly!
%
sq.delay(100e-3);
sq.find(img_ch).set(0);
sq.find('87 repump').set(0);
sq.find('MOT bias').set(0); %ttl on imaging coil

end