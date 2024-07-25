
FigNum = 10;
XLABEL = 'Total Power (mW)';
YLABEL = 'F = 2 Pop';

x_sorted = r.data.Param(1:numel(r.data.R(:,1)));
y_sorted = r.data.R(:,1);

y_sorted(3) = NaN;


[unique_x, ~, idx_uniq] = unique(x_sorted);
averaged_y = accumarray(idx_uniq, y_sorted, [], @mean);


[fitresult, ~] = createFit(unique_x, averaged_y);

x2 = linspace(unique_x(1),unique_x(end),numel(unique_x)*5);
y2 = fitresult(x2);


tau_pi = fitresult.g + (atan(fitresult.b/fitresult.d) + 0*pi - fitresult.c)/fitresult.b;


figure(FigNum);clf
scatter(unique_x, averaged_y); % Example plot
hold on
plot(x2,y2)
xlabel(XLABEL,'Interpreter','latex','FontSize',16)
ylabel(YLABEL,'Interpreter','latex','FontSize',16)

annotation('textbox', [0.25, 0.1, 0.1, 0.1], 'String', sprintf('Tau_pi: %g', tau_pi), ...
    'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');


clear x2 y2 XLABEL YLABEL FigNum x2  x_sorted y_sorted unique_x averaged_y idx_uniq


function [fitresult, gof] = createFit(unique_x, averaged_y)

% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( unique_x, averaged_y );

% Set up fittype and options.
ft = fittype( 'a*sin(b*(x-g)+c)*exp(-d*(x-g))+f', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.859947605330542 0.823835382277067 0.405605888154915 0.0686301996375189 0.610717217585178 0.99098039626571];
opts.Lower = [0 0 0 0 0 -Inf];
opts.Upper = [1 2 Inf Inf 1 Inf];
% opts.Lower = [0 0 -3.2 0 0 0];
% opts.Upper = [1 100 3.2 1000 1 1000];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end
