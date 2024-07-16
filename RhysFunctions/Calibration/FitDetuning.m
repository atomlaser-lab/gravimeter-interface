
% Inputs
x = r.data.Param(1:numel(r.data.R(:,1)));
y = r.data.R(:,1).';
% y = y_av;
x_label = 'Two Photon Detuning (kHz)';
y_label = 'F = 2 Population';
FigNum = 10;

% y(1) = 0;
y(17) = NaN;


[fitresult, ~] = createFit(x, y, 80, max(y), 50);
x2 = linspace(x(1),x(end),numel(x)*10);


figure(FigNum);clf
scatter(x,y)
hold on
plot(x2,fitresult(x2))

annotation('textbox', [0.25, 0.1, 0.1, 0.1], 'String', sprintf('Centre: %g, Width: %g, Amp %g', fitresult.b1, fitresult.c1,fitresult.a1), ...
    'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');
xlabel(x_label,'Interpreter','latex','FontSize',14)
ylabel(y_label,'Interpreter','latex','FontSize',14)

clear x y x_label ylabel FigNum x2 



function [fitresult, gof] = createFit(x, y, CentreGuess, AmpGuess, WidthGuess)

stdGuess = WidthGuess/2;

[xData, yData] = prepareCurveData( x, y );
% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [AmpGuess CentreGuess stdGuess];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end