function makeGaussianPulse(dds,varargin)

%% Set up variables and parse inputs
t0 = 5e-3;
width = 250e-6;
dt = 5e-6;
power = 0.5;
power1 = [];
power2 = [];
df = 153.5e-3;

if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 't0'
                t0 = v;
            case 'dt'
                dt = v;
            case 'width'
                width = v;
            case 'power'
                power = v;
            case 'power1'
                power1 = v;
            case 'power2'
                power2 = v;
            case 'df'
                df = v;
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end

%% Calculate intermediate values
fwhm = width/(2*sqrt(log(2)));

if isempty(power2) && isempty(power1)
    power2 = power;
end

if isempty(power1)
    power1 = power;
end

%% Create vectors
min_t = t0 - 5*width;
max_t = t0 + 5*width;
t = (min_t:dt:max_t)';

P = power1*exp(-(t - t0).^2/fwhm.^2);
P(:,2) = power2*exp(-(t - t0).^2/fwhm.^2);

freq(:,1) = dds(1).DEFAULT_FREQ;
freq(:,2) = freq(:,1) + df;

ph = zeros(numel(t),2);

%% Populate DDS values
for nn = 1:numel(dds)
    dds(nn).after(t,freq(:,nn),P(:,nn),ph(:,nn));
end


end