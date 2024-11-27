function makeNDImagingSequence(sq,varargin)

%
% Define default parameters
%
pulseTime = 30e-6;
camTime = 5e-6; 
cycleTime = 5e-3; 
imgAmplitude = 0.1;
num_images = 1;
pulse_delay = 0;
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
            case 'pulse time'
                pulseTime = v;
            case 'cycle time'
                cycleTime = v;
            case 'cam time'
                camTime = v;
            case 'imaging freq'
                imgFreq = v;
            case 'imaging amplitude'
                imgAmplitude = v;
            case 'num_images'
                num_images = v;
            case 'pulse delay'
                pulse_delay = v;
            otherwise
                error('Unsupported option %s',p);
        end
    end
end
%
% Exit if no images are requested
%
if num_images == 0
    return
end
%
% Preamble - set the imaging frequency and power
%
sq.find('ND imag amp').before(10e-3,imgAmplitude); % Voltage sent to VVA

%
% Imaging beam and camera trigger for image with atoms
%
for nn = 1:num_images
    sq.find('ND imag').after(pulse_delay,1).after(pulseTime,0);     
    sq.find('ND cam trig').set(1).after(camTime,0);                 
    sq.delay(cycleTime);
end
% sq.anchor(sq.latest);   %Re-anchor the sequence to the latest value


end