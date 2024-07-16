
x = r.data.Param(1:numel(r.data.R(:,1)));
y = r.data.R(:,1).';
% x = r.data.Param(1:8);
% y = r.data.R(1:8,1).';

FigNum = 10;
x_label = 'T (ms)';
y_label = 'F = 2 Population';

PeriodEstimate = 50;

[fitresult, gof] = FitSin(x, y,PeriodEstimate,max(y));
x2 = linspace(x(1),x(end),numel(x)*10);


figure(FigNum);clf
scatter(x,y)
hold on
plot(x2,fitresult(x2))
xlabel(x_label,'Interpreter','latex','FontSize',14)
ylabel(y_label,'Interpreter','latex','FontSize',14)
% 
% annotation('textbox', [0.25, 0.1, 0.1, 0.1], 'String', sprintf('Freq: %g, Amp: %g, phi_{off}: %g:, offset: %g', 2*pi/fitresult.b, fitresult.a,fitresult.c,fitresult.d), ...
%     'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');

topaste(FigNum)


function [fitresult, gof] = FitSin(x, y,PeriodEstimate,AmpEstimate)

FrequencyEstimate = 2*pi/PeriodEstimate;

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
% ft = fittype( 'a*sin(b*x+c)+d', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a*sin(x+c)+d', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [AmpEstimate pi/2 0.5];
opts.Lower = [0 -inf -inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end
