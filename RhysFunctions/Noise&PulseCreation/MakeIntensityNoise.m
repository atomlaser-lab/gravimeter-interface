function I_factor = MakeIntensityNoise(varargin)
% When atoms are in free-fall on a moving platform, transverse
% accelerations result in the atoms moving relative to the lasers used.
% This function calculates how the what fraction of the peak intensity the
% atoms will experience for a Gaussian beam.

% The output is a vector of length equal to the number of DDS updates used
% for the specified pulse. This vector should be multiplied to the
% amplitude values of the DDS


% % % Constants
% use beam radius of the long drop
w_0 = 12.5e-3; %("~2 cm FWHM" section 4.2 of Hardman)

% % % Default Values
type = 'acceleration';
amp = 0;
NumPulseWidths = 5;
width = 30e-6;
dt = 1e-6;
t0 = 10e-3;
T = 1e-3;
Tasym = 0;
% a_trans = transverse acceleration
% tof = time since drop that the pulse occurs at
% width = Total pulse duration used by DDS
% dt = dds discritisation of pulse shape


if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 'type'
                type = v;
            case 'amp'
                amp = v;
            case 'width'
                width = v;
            case 'dt'
                dt = v;
            case 't0'
                t0 = v;
            case 't'
                T = v;
            case 'tasym'
                Tasym = v;
            case 'numpulses'
                NumPulses = v;
            case 'numpulsewidths'
                NumPulseWidths = v;
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end



% % % % Create time vector of length equal to the pulse vector
tc = [t0;t0+T;t0+2*T];
for ii =1:length(tc)
    t(:,ii) = (tc(ii)-NumPulseWidths*width:dt:tc(ii)+NumPulseWidths*width);
end
t = t(:);


% % % % Create intensitity noise profile
if strcmpi(type,'acceleration')
    r = 0.5*amp*t.^2;
    I_factor = exp(-2*r.^2/w_0^2);
elseif strcmpi(type,'white')
    I_factor = amp*normrnd(0,1,length(t),1);
elseif strcmpi(type,'white') == 0 && strcmpi(type,'acceleration')
    warning('Type must be "white" or "acceleration"')
    return
end




end