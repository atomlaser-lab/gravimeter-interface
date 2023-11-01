function [NKT_SPD NKT_k_eff] = FrequencyCheck(varargin)
% % % % wavemeter reading:
% Previously Measured Values
Wavemeter_CarrierFreq = 3.842421200000000e14;
Wavemeter_SidebandFreq = 3.842489430000000e+14;

% Set Values
NKT_CarrierWavelength = 1560.4530e-9;
F_Mod_Set = 3.407341305e9;

% Previously Measured Two-Photon Detuning (TPD)
MicrowaveFrequency = 3.417182030452000e09;
TPD_Estimate = 2*MicrowaveFrequency;

%% Constants
WavemeterOffset = -1.65e9;

%Rb87 D2 Transition |F = 2> -> |F' = 0>
F2_FP0 = 384.230484468562e12 - 2.56300597908910934e9 - 302.073888e6;
HFS = 6.83468261090429090e9;


%% Default Inputs
Sideband = 1; %plus or minus integer

% Beam Parameters
w0 = 3e-3;
P1 = 10e-3;
P2 = 10e-3;

% Pulse Parameters
tau = 10e-6;


%% Input Lookup
if mod(numel(varargin),2) ~= 0
    error('Arguments must appear as name/value pairs!');
else
    for nn = 1:2:numel(varargin)
        v = varargin{nn+1};
        switch lower(varargin{nn})
            case 'carrierwavelength'
                NKT_CarrierWavelength = v;
            case 'f_mod'
                F_Mod_Set = v;
            case 'sideband'
                Sideband = v;
            otherwise
                error('Option %s not supported',varargin{nn});
        end
    end
end

%% Calculate Using Measured Frequnecies
% True Frequencies
CarrierFreq = Wavemeter_CarrierFreq + WavemeterOffset;
SidebandFreq = Wavemeter_SidebandFreq + WavemeterOffset;
TPD = abs(SidebandFreq - CarrierFreq);

% Single-Photon Detuning
SPD1 = F2_FP0 - CarrierFreq;
SPD2 = F2_FP0 - SidebandFreq;
SPD = sign(SPD1)*min(abs([SPD1,SPD2]));

% k-vector
k1 = 2*pi*CarrierFreq/const.c;
k2 = 2*pi*SidebandFreq/const.c;
k_eff = abs(k1 - k2);


% % % Calculate RabiFrequency
I1 = P1/(pi*w0^2);
I2 = P2/(pi*w0^2);
d = 1.731e-29;
Coeff = d^2/(const.hbar^2*const.eps0*const.c);

Rabi = abs(Coeff*(sqrt(I1)*sqrt(I2)/SPD));
PulseArea = Rabi*tau;
AproxAOMOffset = (TPD_Estimate - TPD)/1e6/4;


sprintf(['Measured Values\n' ...
    'Single Photon Detuning: %.3g GHz\n' ...
    'Rabi Frequency: %.3g MHz\n' ...
    'k_eff: %.7g\n' ...
    'Pulse Area: %.3g\n' ...
    'AOM Offset: %.5g MHz'], ...
    SPD/1e9,Rabi/1e6,k_eff,PulseArea,AproxAOMOffset)


%% Caclulate Using values stated in NKT 
% Sideband Wavelength
NKT_CarrierFreq = const.c/(NKT_CarrierWavelength);
NKT_SidebandFreq = NKT_CarrierFreq + Sideband*F_Mod_Set;
NKT_SidebandWavelength = const.c/NKT_SidebandFreq;

% Frequency Double
NKT_W1 = NKT_CarrierWavelength/2;
NKT_W2 = NKT_SidebandWavelength/2;

NKT_F1 = const.c/NKT_W1;
NKT_F2 = const.c/NKT_W2;

% Single-Photon Detuning
NKT_SPD1 = F2_FP0 - NKT_F1;
NKT_SPD2 = F2_FP0 - NKT_F2;
NKT_SPD = sign(NKT_SPD1)*min(abs([NKT_SPD1,NKT_SPD2]));

% k-vectors
NKT_k1 = 2*pi*NKT_F1/const.c;
NKT_k2 = 2*pi*NKT_F2/const.c;
NKT_k_eff = abs(NKT_k1-NKT_k2);

NKT_Rabi = abs(Coeff*(sqrt(I1)*sqrt(I2)/NKT_SPD));
NKT_PulseArea = NKT_Rabi*tau;
NKT_AproxAOMOffset = (TPD_Estimate - 2*F_Mod_Set)/1e6/4;
 




%% Outputs
% % % NKT Values 
sprintf(['NKT Values\n' ...
    'Single Photon Detuning: %.3g GHz\n' ...
    'Rabi Frequency: %.3g MHz\n' ...
    'k_eff: %.7g\n' ...
    'Pulse Area: %.3g\n' ...
    'AOM Offset: %.5g MHz'], ...
    NKT_SPD/1e9,NKT_Rabi/1e6,NKT_k_eff,NKT_PulseArea,NKT_AproxAOMOffset)

end