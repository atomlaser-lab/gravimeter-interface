% % % Inputs

FigNum = 10;
% directoryPath = 'E:\R\NoIntensityRamp\DisplacementScan\T scan\T = 2 ms';
% directoryPath = 'E:\R\NoIntensityRamp\DisplacementScan\T scan\T = 20 ms V2';
% legend('T = 2', 'T = 20')

% directoryPath = 'E:\R\NoIntensityRamp\DisplacementScan\t0 scan\t0 = 0 ms';
% directoryPath = 'E:\R\NoIntensityRamp\DisplacementScan\t0 scan\t0 = 30 ms';
% legend('t0 = 0', 't0 = 30')

% directoryPath = 'E:\R\NoIntensityRamp\DisplacementScan\w0 scan\w0  = 1 mm';
% directoryPath = 'E:\R\NoIntensityRamp\DisplacementScan\w0 scan\w0 = 15 mm';
% legend('w0 = 1', 'w0 = 15')

% directoryPath = 'E:\R\NoIntensityRamp\accelScan\T = 2.5 ms';
% directoryPath = 'E:\R\NoIntensityRamp\accelScan\T = 30';
% legend('T = 2.5', 'T = 30')

% directoryPath = 'E:\R\Accelerometer Noise\Bias\T = 10,w0 = 10, t0 = 7.5\T = 10,w0 = 10, t0 = 7.5, Bias = 100e-1';
% directoryPath = 'E:\R\Accelerometer Noise\Bias\T = 10,w0 = 10, t0 = 7.5\T = 10,w0 = 10, t0 = 7.5, Bias = 200e-6';

data = importfile(directoryPath);


ClearFigOnOff = 0;

Title = 'No Ramp';
SubTitle = 't0 Comparison';
XLabel = 'Normalised Displacement';
XLabel = 'acceleration (m/s/s)';

% % % load data
% grab data from file

XX = data.param2(1:end-1);
y = data.raw.R(:,:,1,1);
x = data.param1;
x2 = linspace(0,200,100);

N_Total = data.raw.N(:,1,1,2) + data.raw.N(:,1,1,1);
Nsum_Total = data.raw.Nsum(:,1,1,2) + data.raw.Nsum(:,1,1,1);

for ii = 1:(size(y,2) -1)
% % % Remove outliers
%     Remove N = 0



    warning('off','all')
    [fitresult, ~] = createFit(x, y(:,ii,1,1));
    warning('on','all')

    ci = confint(fitresult);
    a(ii) = fitresult.a;
    b(ii) = fitresult.b;
    c(ii) = fitresult.c;

    a_uncertainty(ii) = (ci(2,1) - ci(1,1)) / 2;
    b_uncertainty(ii) = (ci(2,2) - ci(1,2)) / 2;
    c_uncertainty(ii) = (ci(2,3) - ci(1,3)) / 2; 

%     y_fit(:,ii) = fitresult(x2);
end

figure(FigNum);
if ClearFigOnOff == 1
    clf
end

subplot(3,1,1); hold on
errorbar(XX,a,a_uncertainty,"-s",'MarkerSize',4)
grid on
grid minor
ylabel('Population Offset','Interpreter','latex','FontSize',14)
% ylim([0 0.5])

subplot(3,1,2); hold on
errorbar(XX,b,b_uncertainty,"-s",'MarkerSize',4)
grid on
grid minor
ylabel('Contrast','Interpreter','latex','FontSize',14)
ylim([0,1])

subplot(3,1,3);hold on
errorbar(XX,c - pi/2,c_uncertainty)

grid on
grid minor
ylabel('Phase Offset','Interpreter','latex','FontSize',14)

xlabel(XLabel,'Interpreter','latex','FontSize',14)
sgtitle({['{\bf\fontsize{14}' Title '}'],SubTitle});

clear a b c a_uncertainty b_uncertainty c_uncertainty x x2 XX y
clear N_Total Nsum_Total FigNum ci ClearFigOnOff data ii



function [fitresult, gof] = createFit(x, y)


[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a+b*sin(x*pi/180+c).^2', 'independent', 'x', 'dependent', 'y' );


opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Lower = [0 0 -3.2];
% opts.StartPoint = [0.817627708322262 0.794831416883453 0.644318130193692];
% opts.Upper = [1 1 3.2];

opts.Lower = [0 0 -3.2];
opts.StartPoint = [0.817627708322262 0.794831416883453 0.644318130193692];
opts.Upper = [1 1 3.2];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end






function data = importfile(directoryPath)
% IMPORTFILE Imports the first variable from the 'data.mat' file in the specified directory
% and assigns it to 'data'.
%   directoryPath: Directory path containing the 'data.mat' file.

    % Construct the full path to the 'data.mat' file
    fileToRead = fullfile(directoryPath, 'data.mat');

    % Check if the file exists
    if isfile(fileToRead)
        % Load the .mat file
        newData = load('-mat', fileToRead);
        
        % Get the variable names in the file
        vars = fieldnames(newData);
        
        % Check if there are any variables in the file
        if ~isempty(vars)
            % Assign the first variable found to 'data'
            data = newData.(vars{1});
            
            % Optionally, you can assign this variable to the base workspace as 'data'
            assignin('base', 'data', data);
        else
            error('The file is empty and contains no variables.');
        end
    else
        error('File ''data.mat'' not found in the specified directory: %s', directoryPath);
    end
end