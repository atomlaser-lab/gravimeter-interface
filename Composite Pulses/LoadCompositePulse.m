clear;clc
% % % Inputs
BS_FileName = '20230509T030947_raman_beamsplitter.hdf5';
M_FileName = '20230511T094824_raman_mirror.hdf5';

SPD = 6.8e9;
w0 = 12.5e-3;
I2OnI1 = 1;

FigNum = 10;

%% Import Beam splitter pulse
BSdata = Loadhdf5File(BS_FileName);
Mdata = Loadhdf5File(M_FileName);

% % % This is the data that is used in the "ConvertPulseToDDS" function


%% Convert Raw Data to Rabi Frequency
% Beam Splitter
BS_time = cumsum(BSdata.durs);
BS_duration = BSdata.durs;

BS_TPD = BSdata.dets;
BS_RabiAmp = sqrt(BSdata.rabis.i.^2 + BSdata.rabis.r.^2);
BS_RabiAngle = atan(BSdata.rabis.i ./ BSdata.rabis.r);


% Mirror
M_time = cumsum(Mdata.durs);
M_duration = Mdata.durs;

M_TPD = Mdata.dets;
M_RabiAmp = sqrt(Mdata.rabis.i.^2 + Mdata.rabis.r.^2);
M_RabiAngle = atan(Mdata.rabis.i ./ Mdata.rabis.r);





%% Convert Rabi Frequnecy to Intensity
[BS_MaxLaserPowerRequired,BS_I1Peak,BS_I2Peak] = CalcLaserPower(BS_RabiAmp,SPD,w0,I2OnI1);
[M_MaxLaserPowerRequired,M_I1Peak,M_I2Peak] = CalcLaserPower(M_RabiAmp,SPD,w0,I2OnI1);

%% Convert TPD into an AOM ramp
% There are two double pass AOMS. Thus, a change in frequency of f in the
% AOM results in a 2f change in the beam's frequency.
% Ideally, the two beams are ramped equally such so that the change in
% diffraction efficiency is the same for both beams (thus the desired
% intensity ratio is maintatined)


BS_AOM1FrequencyChange = BS_TPD/4;
BS_AOM2FrequencyChange = -BS_TPD/4;

M_AOM1FrequencyChange = M_TPD/4;
M_AOM2FrequencyChange = -M_TPD/4;


%% Calculate laser phase
BS_Phase = atan2(BSdata.rabis.i,BSdata.rabis.r);
phi2 = atan(BSdata.rabis.i./BSdata.rabis.r);
M_Phase = atan2(BSdata.rabis.i,BSdata.rabis.r);


%%
figure(FigNum+1);clf
subplot(3,1,1)
plot(BS_time*1e6,BS_I2Peak,'r')
hold on
plot(BS_time*1e6,BS_I1Peak,'b')
legend('Carrier', 'Sideband')
ylabel('Beam Power (W)','Interpreter','latex')
xlabel('Time ($\mu s$)', 'Interpreter','latex')

subplot(3,1,2)
plot(BS_time*1e6,BS_AOM1FrequencyChange*1e-6,'b')
hold on
plot(BS_time*1e6,BS_AOM2FrequencyChange*1e-6,'r')
legend('Carrier', 'Sideband')
ylabel('AOM frequency (MHz)','Interpreter','latex')
xlabel('Time ($\mu s$)', 'Interpreter','latex')


subplot(3,1,3)
plot(BS_time*1e6,BS_Phase,'b')
hold on
plot(BS_time*1e6,zeros(size(BS_Phase)),'r')
legend('Carrier', 'Sideband')
ylabel('AOM Phase (rad)','Interpreter','latex')
xlabel('Time ($\mu s$)', 'Interpreter','latex')

%% Plot Pulse
figure(FigNum);clf
subplot(3,2,1)
plot(BS_time*1e6,BS_RabiAmp*1e-3/(2*pi))
ylabel('$|\Omega_R|$ / $2\pi$ (kHz)','Interpreter','latex')
title('Beam Splitter')

grid on
subplot(3,2,3)
plot(BS_time*1e6,BS_RabiAngle)
ylabel('$\Omega _R$ Angle','Interpreter','latex')
grid on

subplot(3,2,5)
plot(BS_time*1e6,BS_TPD*1e-3/(2*pi))
grid on
ylabel('$\delta$ / $2\pi$ (kHz)','Interpreter','latex')
xlabel('Time ($\mu s$)', 'Interpreter','latex')


subplot(3,2,2)
plot(M_time*1e6,M_RabiAmp*1e-3/(2*pi))
ylabel('$|\Omega_R|$ / $2\pi$ (kHz)','Interpreter','latex')
title('Mirror')

grid on
subplot(3,2,4)
plot(M_time*1e6,M_RabiAngle)
ylabel('$\Omega _R$ Angle','Interpreter','latex')
grid on

subplot(3,2,6)
plot(M_time*1e6,M_TPD*1e-3/(2*pi))
grid on
ylabel('$\delta$ / $2\pi$ (kHz)','Interpreter','latex')
xlabel('Time ($\mu s$)', 'Interpreter','latex')



























%% functions
function data = Loadhdf5File(FileName)
% % % % % h5disp(FileName) %This prints the data info
fileInfo = h5info(FileName);
LevelName = fileInfo.Groups.Groups(5).Name;
DataSetName = fileInfo.Groups.Groups(5).Datasets.Name;
PathName = append(LevelName,'/',DataSetName);
data = h5read(FileName,PathName);
end



function [MaxLaserPowerRequired,I1Peak,I2Peak] = CalcLaserPower(RabiMax,SinglePhotonDetuning,w0,I2onI1)
% RabiMax = d^2/(h^2*e*c) * (I/Delta)
TransitionDipoleMoment = 1.731e-29; % C m
Coeff = (TransitionDipoleMoment^2)/(const.hbar^2*const.eps0*const.c*SinglePhotonDetuning);

I2Peak = sqrt((RabiMax.^2/Coeff^2)*I2onI1);
I1Peak = I2Peak/I2onI1;

P1Peak = I2Peak*(pi*w0^2);
P2Peak = I1Peak*(pi*w0^2);

MaxLaserPowerRequired = max(P1Peak+P2Peak);
end