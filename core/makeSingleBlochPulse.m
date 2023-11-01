function makeSingleBlochPulse(dds, varargin)

%% Set up variables and parse inputs
f = 384.224e12;
k = 2*pi*f/const.c;
t0 = 10e-3;
rampTime = 10e-3;
power = 0.05;
dt = 1e-6;
mirrorSwitch = 1;

% Parse input parameters
if mod(numel(varargin), 2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 't0'
                t0 = v;
            case 'ramptime'
                rampTime = v;
            case 'power'
                power = v;
            % ... other options can go here
            otherwise
                error('Option %s not supported', varargin{nn});
        end
    end
end

%% Time Vector
t = (0:dt:rampTime)';

%% Calculate intermediate values
recoil = const.hbar * k^2 / (2 * const.mRb * 2 * pi);

%% Initialize powers, phases, and frequencies
[P, ph, freq] = deal(zeros(numel(t), 2));

%% Linear Frequency Ramp for Bloch Oscillation
deltaFreq = recoil / rampTime;  % Frequency change per unit time
freqRamp = linspace(0, deltaFreq, numel(t));

%% Set powers, phases, and frequencies
P(:, 1) = power;
P(:, 2) = power;

freq(:, 1) = DDSChannel.DEFAULT_FREQ + freqRamp;
freq(:, 2) = DDSChannel.DEFAULT_FREQ - freqRamp;

ph(:, :) = 0;

%% Populate DDS values
for nn = 1:numel(dds)
    dds(nn).after(t + t0, freq(:, nn), P(:, nn), ph(:, nn));
end

end
