w0 = 3e-3;
(pi^2*w0^2*const.hbar^2*const.eps0*const.c*1e10)/(1e-6*1.731e-29^2)




% 1.731e-29




%%
% % % Inputs
PulseChain = [0.5,1,0.5];
T = 10e-3;
t0 = 5e-3;

% % % Load Pulse Info
% Load raw data
load RamanPulses_raw.mat

% Calculate DDS instructions based on max power/single-photon detuning/beam size/intensity ratio
[dt1 BS_P1 BS_P2 BS_F1 BS_F2 BS_phi1 BS_phi2] = ConvertPulseToDDS('BS',P1Max,P2Max,SPD,w0,I2onI1);
[dt2 M_P1 M_P2 M_F1 M_F2 M_phi1 M_phi2] = ConvertPulseToDDS('M',P1Max,P2Max,SPD,w0,I2onI1);



% % % Create time vector
numPulses = numel(PulseChain);
MirrorPulseNumInstructions = 400;
TimeVectorLength = PulseChain*MirrorPulseNumInstructions;

% predefine vectors
t = zeros(MirrorPulseNumInstructions,numPulses);
P1 = zeros(MirrorPulseNumInstructions,numPulses);
P2 = zeros(MirrorPulseNumInstructions,numPulses);
F1 = zeros(MirrorPulseNumInstructions,numPulses);
F2 = zeros(MirrorPulseNumInstructions,numPulses);
phi1 = zeros(MirrorPulseNumInstructions,numPulses);
phi2 = zeros(MirrorPulseNumInstructions,numPulses);

% define pulse durations for each pulse type
tau1 = sum(BSdata.durs);
tPulse1 = (-tau1/2:dt1:tau1/2)';

tau2 = sum(Mdata.durs);
tPulse2 = (-tau2/2:dt2:tau2/2);


for nn = 1:numPulses
    if PulseChain(nn) == 0.5
        t(1:length(tPulse1),nn) = tPulse1;
        t(length(tPulse1)+1:end,nn) = NaN;
    elseif PulseChain(nn) == 1
        t(:,nn) = tPulse2;
    else
        warning('Composite pulse must be a pi/2 or pi pulse')
    end
    
    t(:,nn) = t(:,nn) + t0 + T*nn;
end

t = t(:)';
t(isnan(t)) = [];

