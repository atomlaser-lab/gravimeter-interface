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
repumpShutterDelay = 2e-3;
liquidCrystalDelay = 12e-3;

% % On/Off Options
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

sq.anchor(timeAtDrop);
sq.find('imaging freq').after(tof1 - 1e-3,imgFreq);
sq.find('Liquid Crystal Repump').after(tof1 - 1e-3,7);
sq.find('Repump Amp').after(tof1 - 1e-3,9);
sq.find('Repump Freq').after(tof1 - 1e-3,RunConversions.repump_freq(0));


% Image atoms in the F = 2 manifold
sq.anchor(timeAtDrop);
sq.find('imaging amp ttl').after(tof1+pulseDelay,1).after(pulseTime,0);
sq.find(camChannel).after(tof1,1).after(camTime,0);

% set second imaging detuning with lots of time for VCO to change
sq.find('imaging freq').set(imgFreq2);

% Blow away F=2 atoms
% tBlow = tof2 - tof1 - 0.1e-3;
if (tof2 - tof1) < 2e-3
    warning('Difference between TOF1 and TOF2 is less than 2 ms, may interfere with F = 2 blow away pulse');
end
tBlow = 2e-3;
blowDelay = 10e-6;
sq.find('3D MOT Amp TTL').after(tof1 + camTime + blowDelay,1).after(tBlow,0);
% sq.find('3D MOT Freq').set(RunConversions.mot_freq(0)).after(tBlow,RunConversions.mot_freq(-70));
% sq.find('3D MOT Amp').set(5).after(tBlow,0);
sq.find('3D MOT Freq').after(tof1 + camTime + blowDelay,RunConversions.mot_freq(0)).after(tBlow,RunConversions.mot_freq(-70));
sq.find('3D MOT Amp').after(tof1 + camTime + blowDelay,5).after(tBlow,0);

sq.delay(RepumpDelay);



% Image atoms in the F = 1 state
sq.anchor(timeAtDrop);
sq.find('Imaging amp ttl').after(tof2+pulseDelay,1).after(pulseTime,0);
sq.find(camChannel).after(tof2,1).after(camTime,0);
sq.find('repump amp ttl').after(tof2 - pulseDelay,0).before(repumpTime,1);
sq.find('Top repump shutter').after(tof2 - pulseDelay - repumpShutterDelay,0).after(repumpTime+repumpShutterDelay,1);
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