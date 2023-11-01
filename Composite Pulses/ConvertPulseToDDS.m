function [dt P1 P2 F1 F2 phi1 phi2] = ConvertPulseToDDS(PulseType,P1Max,P2Max,SPD,w0,I2onI1)
% Q-CTRL has given us composite pulses that are characterised by the real
% and imaginary component of the pulse's rabi-frequency, the two-photon 
% detuning, and the duration of each step

% This function does a few things
%   1) it converts the real/imaginary components into an amplitude/angle
%   2) it converts the quantities into an intensity/frequency/phase
%   3) it converts the I/f/phi into AOM values

% Script inputs:
%   PMaxi is the maximum available power in beam i
%       i = 1 is the carrier
%   SPD = single photon detuning in GHz
%   w0 = 1/e^2 beam radius in m
%   I2onI1 is the intensity ratio (sideband divided by carrier)
%   PulseData is the loaded pulse data provided by   Q-CTRL



%% Load the pulse data
load RamanPulses_raw.mat
if strcmpi(PulseType,'BS') == 1
    PulseData = BSdata;
elseif strcmpi(PulseType,'M') == 1
    PulseData = BSdata;
end

%% 1: Convert components into amplitude/angle
Time = cumsum(PulseData.durs);
Duration = PulseData.durs;

TPD = PulseData.dets;
RabiAmp = sqrt(PulseData.rabis.i.^2 + PulseData.rabis.r.^2);
RabiAngle = atan(PulseData.rabis.i ./ PulseData.rabis.r);


%% 2: Convert to intensity/phase
% % % Intensity
TransitionDipoleMoment = 1.731e-29; % C m
Coeff = (TransitionDipoleMoment^2)/(const.hbar^2*const.eps0*const.c*SPD);

I2Peak = sqrt((RabiAmp.^2/Coeff^2)*I2onI1);
I1Peak = I2Peak/I2onI1;

P1Peak = I2Peak*(pi*w0^2);
P2Peak = I1Peak*(pi*w0^2);


if PMax1 <= max(P1Peak) || PMax2 <= max(P2Peak)
    warning('More power is required than there is avaialbe for these pulses')
end


% % % Phase
Phase = atan2(PulseData.rabis.i,PulseData.rabis.r);

%% 3: Convert I/f/phi into AOM values
% Note that subscript i indicates AOM i


% % % Pulse timing
% Current pulse builder requires that instructions have the same length

if sum(diff(Duration) == 0) ~= length(Duration) -1 
    dt = Duration(1);
else
    warning('Instructions vary length')
end

% % % Freq
% There are two double pass AOMS. Thus, a change in frequency of f in the
% AOM results in a 2f change in the beam's frequency.
% Ideally, the two beams are ramped equally so that the change in
% diffraction efficiency is the same for both beams (thus the desired
% intensity ratio is maintatined). Hence:
F1 = TPD/4;
F2 = -TPD/4;

% % % Intensity
% makePulseSequence requires a percentage of total power
P1 = P1Peak/PMax1;
P2 = P2Peak/PMax2;

% % % Phase
% Hold one laser constant, and ramp the other
phi1 = zeros(size(Phase));
phi2 = Phase;

end