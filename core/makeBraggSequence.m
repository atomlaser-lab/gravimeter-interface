function makeBraggSequence(dds,varargin)

%% Set up variables and parse inputs
f = 384.224e12;
k = 2*pi*f/const.c;
t0 = 10e-3;
width = 30e-6;
T = 1e-3;
Tasym = 0;
appliedPhase = 0;
power = 0.05*[1,2,1];
power1 = [];power2 = [];
chirp = 2*k*9.795/(2*pi);
order = 1;
start_order = 0;
mirrorSwitch = 1;

if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 't0'
                t0 = v;
            case 't'
                T = v;
            case 'dt'
                dt = v;
            case 'tasym'
                Tasym = v;
            case 'width'
                width = v;
            case {'appliedphase','phase'}
                appliedPhase = v;
            case 'power'
                power = v;
                if any(power < 0)
                    error('Power needs a value between 0 and 1.');
                elseif power > 1
                    error('Power needs a value between 0 and 1.');
                else
                    power = v;
                end
            case 'chirp'
                chirp = v;
            case 'f'
                f = v;
                k = 2*pi*f/const.c;
            case 'k'
                k = v;
            case 'power1'
                power1 = v;
            case 'power2'
                power2 = v;
            case 'order'
                if v == 0
                    error('Bragg order should be non-zero!');
                elseif round(v) ~= v
                    error('Bragg order must be an integer!');
                else
                    order = v;
                end
            case 'start_order'
                 if round(v) ~= v
                     error('Bragg start order must be an integer!');
                else
                    start_order = v;
                end
            case 'mirror'
                switch lower(v)
                    case 'isolated'
                        mirrorSwitch = 1;
                    case {'mot','rigid'}
                        mirrorSwitch = -1;
                    otherwise
                        error('''Mirror'' option can only be either ''isolated'' or ''mot''');
                end
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end

%% Conditions on the time step and the Bragg order
if width<50e-6
    dt=1e-6;
else
    intermediatewidth=width*10^6;
    dt = ceil(intermediatewidth/50)*10^-6;
end
     
%% Calculate intermediate values
recoil = const.hbar*k^2/(2*const.mRb*2*pi);
numPulses = numel(power);
fwhm = width/(2*sqrt(log(2)));
detuning=start_order*const.hbar*k^2/const.mRb;

if isempty(power1)
    power1 = power;
end
if isempty(power2)
    power2 = power;
end

if numel(appliedPhase) == 0
    tmp = zeros(1,numPulses);
    tmp(end) = appliedPhase;
    appliedPhase = tmp;
end

%% Create vectors
tPulse = (-5*width:dt:5*width)';
t = repmat(tPulse,1,numPulses);
for  nn = 1:numPulses
    t(:,nn) = t(:,nn) + t0 + (nn-1)*T + max((nn-2),0)*Tasym;
end
t = t(:);
%
% Set powers, phases, and frequencies
%
[P,ph,freq] = deal(zeros(numel(t),2));
for nn = 1:numPulses
    tc = t0 + (nn-1)*T + max((nn-2),0)*Tasym;
    idx = (t - t0) > (nn-1-0.5)*T;
    %
    % Set powers
    %
    P(:,1) = P(:,1) + power1(nn)*exp(-(t - tc).^2/fwhm.^2);
    P(:,2) = P(:,2) + power2(nn)*exp(-(t - tc).^2/fwhm.^2);
    %
    % Set phases
    %
    ph(idx,2) = appliedPhase(nn);
    %
    % Set frequencies.  Need channel 1 frequency to be higher than channel
    % 2 frequency so that we use the lattice formed from retroreflecting
    % from the vibrationally isolated mirror
    %
%     freq(idx,1) = DDSChannel.DEFAULT_FREQ + mirrorSwitch*0.25/1e6*(chirp*tc + order*4*recoil);
%     freq(idx,2) = DDSChannel.DEFAULT_FREQ - mirrorSwitch*0.25/1e6*(chirp*tc + order*4*recoil);
    
    freq(idx,1) = DDSChannel.DEFAULT_FREQ + mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order)*4*recoil);
    freq(idx,2) = DDSChannel.DEFAULT_FREQ - mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order)*4*recoil); 
%    freq(idx,1) = 20;
%     freq(idx,2) = 20;


%     freq(idx,1) = DDSChannel.DEFAULT_FREQ + mirrorSwitch*0.25/1e6*(chirp*tc + (order+detuning/recoil)*4*recoil);
%     freq(idx,2) = DDSChannel.DEFAULT_FREQ - mirrorSwitch*0.25/1e6*(chirp*tc + (order+detuning/recoil)*4*recoil);

%     freq(idx,1) = DDSChannel.DEFAULT_FREQ + mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order+detuning/recoil)*4*recoil);
%     freq(idx,2) = DDSChannel.DEFAULT_FREQ - mirrorSwitch*0.25/1e6*(chirp*tc + (2*start_order+order+detuning/recoil)*4*recoil);

end

freq(freq == 0) = DDSChannel.DEFAULT_FREQ;

%% Populate DDS values
for nn = 1:numel(dds)
    dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
end


end