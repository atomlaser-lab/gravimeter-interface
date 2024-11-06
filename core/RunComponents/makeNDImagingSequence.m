function makeNDImagingSequence(sq,varargin)

%
% Define default parameters
%
pulseTime = 30e-6;
camTime = 5e-6; 
cycleTime = 5e-3; 
imgFreq = 8.5;
imgAmplitude = 1e-3; % 1mW of power 
imgPower   = imgAmplitude; %% addition
imgVoltage = TrapPtoV_NDI(imgPower);
num_images = 1;
species = 85;
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
                % imgAmplitude = v;
                imgPower = v;
            case 'imaging voltage' % ADDED
                imgVoltage = v;
            case 'num_images'
                num_images = v;
            case 'species'
                species = v;
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
% Preamble - set the imaging frequency
%
if species == 87
    sq.find('87 imag freq').set(imgFreq);
    sq.find('87 imag amp').set(TrapPtoV('nd',imgAmplitude));
    ch = sq.find('87 imag');
elseif species == 85
    sq.find('85 imag freq').set(imgFreq);
    sq.find('85 imag amp').set(TrapPtoV_NDI(imgPower)); % power in W
    % sq.find('85 imag amp').set(imgVoltage); % Voltage sent to VVA
    ch = sq.find('85 imag');
end


%
% Imaging beam and camera trigger for image with atoms
%
for nn = 1:num_images
    ch.after(pulse_delay,1).after(pulseTime,0);     %Turn on after TOF, then turn off after pulse time
    sq.find('ND cam trig').set(1).after(camTime,0);   %Turn on after TOF, then turn off after camera time
    sq.delay(cycleTime);
end
% sq.anchor(sq.latest);   %Re-anchor the sequence to the latest value


end