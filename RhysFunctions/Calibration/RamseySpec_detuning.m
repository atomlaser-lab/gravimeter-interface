clear x y x2 

x = r.data.Param(1:numel(r.data.R(:,1)));
y = r.data.R(:,1).';
FigNum = 5;
x_label = 'Pulse Separation Time (ms)';
y_label = 'F = 2 Population';

PeriodEstimate = 6;

[fitresult, gof] = FitSin(x, y,PeriodEstimate,max(y));
x2 = linspace(x(1),x(end),numel(x)*10);


figure(FigNum);clf
scatter(x,y)
hold on
plot(x2,fitresult(x2))
xlabel(x_label,'Interpreter','latex','FontSize',14)
ylabel(y_label,'Interpreter','latex','FontSize',14)

annotation('textbox', [0.25, 0.1, 0.1, 0.1], 'String', sprintf('Freq: %g, Amp: %g', 2*pi/fitresult.b, fitresult.a), ...
    'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');


topaste(5)


function [fitresult, gof] = FitSin(x, y,PeriodEstimate,AmpEstimate)

FrequencyEstimate = 2*pi/PeriodEstimate;

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a*sin(b*x+c)+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [AmpEstimate FrequencyEstimate pi/2 0.5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end


