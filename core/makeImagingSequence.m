function makeImagingSequence(sq,varargin)

imgType = 'in-trap';
pulseTime = [];
repumpTime = 100e-6;
repumpDelay = 00e-6;
fibreSwitchDelay = 20e-3;
camTime = 100e-6;
pulseDelay = 0;
cycleTime = 100e-3;
repumpFreq = 4.3;
imgFreq = 8.5;
manifold = 1;
includeDarkImage = false;
if mod(numel(varargin),2) ~= 0
    error('Input arguments must be in name/value pairs');
else
    for nn = 1:2:numel(varargin)
        p = lower(varargin{nn});
        v = varargin{nn+1};
        switch p
            case 'tof'
                tof = v;
            case 'type'
                imgType = v;
            case 'pulse time'
                pulseTime = v;
            case 'repump time'
                repumpTime = v;
            case 'repump delay'
                repumpDelay = v;
            case 'pulse delay'
                pulseDelay = v;
            case 'cycle time'
                cycleTime = v;
            case 'cam time'
                camTime = v;
            case 'repump freq'
                repumpFreq = v;
            case 'imaging freq'
                imgFreq = v;
            case 'fibre switch delay'
                fibreSwitchDelay = v;
            case 'manifold'
                manifold = v;
            case 'includedarkimage'
                includeDarkImage = v;
            otherwise
                error('Unsupported option %s',p);
        end
    end
end

switch lower(imgType)
    case {'in trap','in-trap','trap','drop 1'}
        camChannel = 'cam trig';
        imgType = 0;
        if isempty(pulseTime)
            pulseTime = 13e-6; % was 17e-6 on the 19/11/23
        end
    case {'drop 2'}
        camChannel = 'drop 1 camera trig';
        imgType = 1;
        if isempty(pulseTime)
            pulseTime = 50e-6;
        end
    otherwise
        error('Unsupported imaging type %s',imgType);
end

%Preamble
sq.find('imaging freq').set(imgFreq);
sq.find('Top repump shutter').before(10e-3,0);

%Repump settings - repump occurs just before imaging
%If manifold is set to image F = 1 state, enable repump. Otherwise,
%disable repumping
if imgType == 0 && manifold == 1
    sq.find('liquid crystal repump').set(7);%was -2.22
    sq.find('repump amp ttl').after(tof-repumpTime-repumpDelay,1);
    sq.find('repump amp ttl').after(repumpTime,0);
    if ~isempty(repumpFreq)
        sq.find('repump freq').after(tof-repumpTime-repumpDelay,repumpFreq);
    end
elseif imgType == 1 && manifold == 1
    sq.find('liquid crystal repump').set(7);
    sq.find('drop repump').after(tof-repumpTime-repumpDelay,1);
    sq.find('drop repump').after(repumpTime,0);
    sq.find('fiber switch repump').after(tof-fibreSwitchDelay,1);   
    if ~isempty(repumpFreq)
        sq.find('drop repump freq').after(tof-repumpTime-repumpDelay,4.3);
    end
end

%Imaging beam and camera trigger for image with atoms
sq.find('Imaging amp ttl').after(tof+pulseDelay,1);
sq.find(camChannel).after(tof,1);
sq.find('imaging amp ttl').after(pulseTime,0);
sq.find(camChannel).after(camTime,0);
sq.anchor(sq.latest);
sq.delay(cycleTime);

%Take image without atoms
sq.find('Imaging amp ttl').after(pulseDelay,1);
sq.find(camChannel).set(1);
sq.find('imaging amp ttl').after(pulseTime,0);
sq.find(camChannel).after(camTime,0);
sq.anchor(sq.latest);
sq.find('fiber switch repump').set(0);

if includeDarkImage
    %Take dark image
    sq.delay(cycleTime);
    sq.find('Imaging amp ttl').after(pulseDelay,0);
    sq.find(camChannel).set(1);
    sq.find('imaging amp ttl').after(pulseTime,0);
    sq.find(camChannel).after(camTime,0);
    sq.anchor(sq.latest);
    sq.find('fiber switch repump').set(0);
end



end