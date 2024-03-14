x = Dipole.W50.Voltage(1:end);
y = Dipole.W50.Power(1:end);


[fitresult, gof] = createFit(x, y);
xtemp = (0:0.2:27.6);

figure(5);clf
subplot(2,1,1)
scatter(y,x)
hold on
plot(y,fitresult(y))


subplot(2,1,2)
scatter(y,fitresult(y) - x)

function [fitresult, gof] = createFit(x, y)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( y, x );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
end