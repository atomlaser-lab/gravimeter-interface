function makeImagingSequence(sq,varargin)

%
% Define default parameters
%
pulseTime = 30e-6;
repumpTime = 100e-6;
repumpDelay = 00e-6;
fibreSwitchDelay = 20e-3;
camTime = 100e-6;
cycleTime = 40e-3;
repumpFreq = 4.58;
repumpAmplitude = 8;
imgFreq = 8.5;
imgAmplitude = 10;
manifold = 1;
take_dark_image = true;
imaging_field = 2.95;
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
            case 'pulse time'
                pulseTime = v;
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
            case 'fibre switch delay'
                fibreSwitchDelay = v;
            case 'manifold'
                manifold = v;
            case 'take dark image'
                take_dark_image = v;
            case 'imaging_field'
                imaging_field = v;
            case 'image type'
                image_type = v;
            otherwise
                error('Unsupported option %s',p);
        end
    end
end


% Set imaging parameters BEFORE you take the image
sq.find('CD0 Fast').set(0); %zero mag field
sq.find('MOT bias coil').set(imaging_field); %turn on the imaging coil (to align the axis of atoms)
sq.find('MOT bias').set(1); %ttl on imaging coil


%
% Preamble - set the imaging frequency
%
sq.find('87 imag freq').set(imgFreq);
sq.find('87 imag amp').set(imgAmplitude);

%Repump settings - repump occurs just before imaging
%If manifold is set to image F = 1 state, enable repump. Otherwise,
%disable repumping
if manifold == 1
    % Set repump amplitude and frequency
    sq.find('87 repump freq').set(repumpFreq);
    sq.find('87 repump amp').set(repumpAmplitude);
    %Turn on the repump TTL and the fiber switch (inverted!)
    sq.find('87 repump').after(tof-repumpTime-repumpDelay,1).after(repumpTime,0);
    sq.find('Repump Switch').after(tof - fibreSwitchDelay,0);
end

if strcmpi(image_type,'horizontal')
    cam_trig = '87 cam trig';
elseif strcmpi(image_type , 'vertical')
    cam_trig = 'vertical cam trig';
elseif strcmpi(image_type,'85')
    cam_trig = 'ND cam trig';
elseif strcmpi(image_type,'MOT')
    cam_trig = '87 cam trig';
else
    warning('incompatible cam trig input')
end

%
% Imaging beam and camera trigger for image with atoms
%
sq.find('87 imag').after(tof,1).after(pulseTime,0); %Turn on after TOF, then turn off after pulse time
sq.find(cam_trig).after(tof,1).after(camTime,0);    %Turn on after TOF, then turn off after camera time
sq.waitFromLatest(cycleTime);                       %Delay
%
% Take image without atoms
%
sq.find('87 imag').set(1).after(pulseTime,0);       %Turn on after TOF, then turn off after pulse time
sq.find(cam_trig).set(1).after(camTime,0);          %Turn on after TOF, then turn off after pulse time
sq.waitFromLatest(cycleTime);                       %Delay
sq.find('Repump Switch').before(50e-6,1);    %Turn off fiber switch

%
% Take a dark image
%
if take_dark_image
    sq.find(cam_trig).set(1).after(camTime,0);   %Turn on after TOF, then turn off after camera time
    sq.anchor(sq.latest);   %Re-anchor the sequence to the latest value
end
%
% Need a last instruction so that the run ends properly!
%
sq.delay(100e-3);
sq.find('87 imag').set(0);
sq.find('87 repump').set(0);
sq.find('MOT bias').set(0); %ttl on imaging coil

end