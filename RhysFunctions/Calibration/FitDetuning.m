
% Inputs
% x = r.data.Param(1:numel(r.data.R(:,1)));
% y = r.data.Rsum(:,1).';

x = r.data.df(1:numel(r.data.N(:,1)));
y = r.data.C1.Nsum./(r.data.C1.Nsum + r.data.C2.Nsum);

y(1:2) = NaN;
y(end) = NaN;

% x_label = 'Two Photon Detuning (kHz)';
% y_label = 'F = 2 Population';

x_label = 'Chirp';
y_label = '|p_0> Pop';
FigNum = 10;

% y(10) = NaN;
% y(end-2:end) = NaN;


[fitresult, ~] = createFit(x, y, 2.511, 0.6, 6e-5);
x2 = linspace(min(x),max(x),numel(x)*10);


figure(FigNum);clf
scatter(x,y,'r','filled')
hold on
plot(x2,fitresult(x2),'r')

annotation('textbox', [0.25, 0.1, 0.1, 0.1], 'String', sprintf('Centre: %g, Width: %g, Amp %g', fitresult.b1, fitresult.c1,fitresult.a1), ...
    'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');
% 
% annotation('textbox', [0.25, 0.8, 0.1, 0.1], 'String', sprintf('Centre: %g, Width: %g, Amp %g', fitresult.b1, fitresult.c1,fitresult.a1), ...
%     'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');


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