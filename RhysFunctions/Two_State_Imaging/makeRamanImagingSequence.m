function makeRamanImagingSequence(sq,varargin)
% % % Default Values
imgType = 'in-trap';
pulseTime = [];
repumpTime = 100e-6;
fibreSwitchDelay = 20e-3;
camTime = 100e-6;
pulseDelay = 0;
repumpFreq = 4.3;
imgFreq = 8.5;
imgFreq2 = imgFreq;
manifold = 1;
includeDarkImage = true;

RepumpDelay = 1e-3;
cycleTime = 100e-3;

% % On/Off Options
BlowAway = 0;
RepumpOnOff = 1;


if mod(numel(varargin),2) ~= 0
    error('Input arguments must be in name/value pairs');
else
    for nn = 1:2:numel(varargin)
        p = lower(varargin{nn});
        v = varargin{nn+1};
        switch p
            case 'tof1'
                tof1 = v;
            case 'tof2'
                tof2 = v;
            case 'type'
                imgType = v;
            case 'pulse time'
                pulseTime = v;
            case 'repump time'
                repumpTime = v;
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
            case 'image freq2'
                imgFreq2 = v;
            case 'fibre switch delay'
                fibreSwitchDelay = v;
            case 'manifold'
                manifold = v;
            case 'includedarkimage'
                includeDarkImage = v;
            case 'repump delay'
                RepumpDelay = v;
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
            pulseTime = 20e-6;
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
timeAtDrop = sq.time;
sq.find('imaging freq').set(imgFreq);
sq.find('Top repump shutter').before(10e-3,0);



% Image atoms in the F = 2 manifold
sq.anchor(timeAtDrop);
sq.find('Imaging amp ttl').after(tof1+pulseDelay,1);
sq.find(camChannel).after(tof1,1);
sq.find('imaging amp ttl').after(pulseTime,0);
sq.find(camChannel).after(camTime,0);

sq.anchor(sq.latest);
sq.delay(RepumpDelay);

% Blow away atoms
if BlowAway == 1
    sq.find('Imaging amp ttl').set(1);
    sq.delay(2.5e-3);
    sq.find('Imaging amp ttl').set(0);
    sq.delay(1e-3);
end

% set second imaging detuning
sq.find('imaging freq').set(imgFreq2);

% Pump atoms in F = 1 into the F = 2 manifold
if RepumpOnOff == 1
    sq.find('liquid crystal repump').before(tof1,7);%was -2.22
    sq.find('repump amp ttl').set(1);
    sq.find('Repump Amp').set(9);
    sq.find('Repump Freq').set(RunConversions.repump_freq(0));
    sq.find('repump amp ttl').after(repumpTime,0);
    sq.find('liquid crystal repump').after(repumpTime,-2.22);
end
sq.anchor(sq.latest);
sq.delay(repumpTime);

% Image atoms in the F = 1 state
sq.anchor(timeAtDrop);
sq.find('Imaging amp ttl').after(tof2+pulseDelay,1);
sq.find(camChannel).after(tof2,1);
sq.find('imaging amp ttl').after(pulseTime,0);
sq.find(camChannel).after(camTime,0);

sq.anchor(sq.latest);
sq.delay(cycleTime);

% Image without atoms
sq.find('Imaging amp ttl').after(pulseDelay,1);
sq.find(camChannel).set(1);
sq.find('imaging amp ttl').after(pulseTime,0);
sq.find(camChannel).after(camTime,0);
sq.anchor(sq.latest);
sq.find('fiber switch repump').set(0);

% Dark Image
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