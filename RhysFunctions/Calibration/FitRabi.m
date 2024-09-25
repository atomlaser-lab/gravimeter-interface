
FigNum = 10;
XLABEL = 'Total Bragg beam(mW)';
% XLABEL = 'Pulse Duration (us)';
YLABEL = '|p_0 + nhk> Pop';
tau_piEstimate = 100;


% x_sorted = r.data.Param(1:numel(r.data.R(:,1)));
% y_sorted = r.data.Rsum(:,1);


x_sorted = r.data.param(1:numel(r.data.C2.N));
y_sorted = r.data.C2.Nsum./(r.data.C1.Nsum + r.data.C2.Nsum);

y_sorted(end) = NaN;

% Average
[unique_x, ~, idx_uniq] = unique(x_sorted);
averaged_y = accumarray(idx_uniq, y_sorted, [], @mean);

% % Use all data
% unique_x =x_sorted;
% averaged_y = y_sorted;

[fitresult, ~] = createFit(unique_x, averaged_y,tau_piEstimate,max(averaged_y));

x2 = linspace(unique_x(1),unique_x(end),numel(unique_x)*50);
y2 = fitresult(x2);


tau_pi = fitresult.g + (atan(fitresult.b/fitresult.d) + 0*pi - fitresult.c)/fitresult.b;
tau_pi = x2(find(y2 == max(y2)));

figure(FigNum);clf
scatter(unique_x, averaged_y); % Example plot
hold on
plot(x2,y2)
xlabel(XLABEL,'Interpreter','latex','FontSize',16)
ylabel(YLABEL,'Interpreter','latex','FontSize',16)
scatter(tau_pi,fitresult(tau_pi),'filled','k','LineWidth',2)
annotation('textbox', [0.25, 0.1, 0.1, 0.1], 'String', sprintf('Tau_pi: %g', tau_pi), ...
    'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');


clear x2 y2 XLABEL YLABEL FigNum x2  x_sorted y_sorted unique_x averaged_y idx_uniq


function [fitresult, gof] = createFit(unique_x, averaged_y,periodEstimate,MaxEstiamte)

% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( unique_x, averaged_y );

% Set up fittype and options.
ft = fittype( 'a*sin(b*(x-g)+c)*exp(-d*(x-g))+f', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [0.859947605330542 0.823835382277067 0.405605888154915 0.0686301996375189 0.610717217585178 0.99098039626571];
% opts.StartPoint = [0.9 pi/20 0.41 0.069 0.6 1];
opts.StartPoint = [MaxEstiamte pi/periodEstimate 5 0.011 0.53 7];

opts.Lower = [0 0 0 0 -Inf 0];
opts.Upper = [1 1000 Inf Inf Inf 1];
% opts.Lower = [0 0 -3.2 0 0 0];
% opts.Upper = [1 100 3.2 1000 1 1000];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end



% 
% FigNum = 10;
% XLABEL = 'Power in each Bragg beam(mW)';
% % XLABEL = 'Pulse Duration (us)';
% YLABEL = 'F = 2 Pop';
% tau_piEstimate = 40;
% 
% x_sorted = r.data.param(1:numel(r.data.C2.N));
% y_sorted = r.data.C2.N./(r.data.C1.N + r.data.C2.N);
% 
% [unique_x, ~, idx_uniq] = unique(x_sorted);
% averaged_y = accumarray(idx_uniq, y_sorted, [], @mean);
% 
% [fitresult, ~] = createFit(unique_x, averaged_y, tau_piEstimate, max(averaged_y));
% 
% x2 = linspace(unique_x(1), unique_x(end), numel(unique_x) * 50);
% y2 = fitresult(x2);
% 
% % Find the peak of the fitted curve
% [max_y, idx_max] = max(y2);
% tau_pi = x2(idx_max);
% 
% figure(FigNum); clf
% scatter(unique_x, averaged_y); % Example plot
% hold on
% plot(x2, y2)
% xlabel(XLABEL, 'Interpreter', 'latex', 'FontSize', 16)
% ylabel(YLABEL, 'Interpreter', 'latex', 'FontSize', 16)
% scatter(tau_pi, max_y, 'filled', 'k', 'LineWidth', 2)
% annotation('textbox', [0.25, 0.1, 0.1, 0.1], 'String', sprintf('Tau_pi: %g', tau_pi), ...
%     'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');
% 
% clear x2 y2 XLABEL YLABEL FigNum x2 x_sorted y_sorted unique_x averaged_y idx_uniq
% 
% function [fitresult, gof] = createFit(unique_x, averaged_y, periodEstimate, MaxEstimate)
% 
% % Fit: 'fit with logistic constraint'.
% [xData, yData] = prepareCurveData(unique_x, averaged_y);
% 
% % Define the logistic function
% ft = fittype(@(a, b, c, d, e, x) 1 ./ (1 + exp(-(a * x + b) / c)) * (d - e) + e, ...
%     'independent', 'x', 'dependent', 'y');
% 
% % Set up fit options
% opts = fitoptions('Method', 'NonlinearLeastSquares');
% opts.Display = 'Off';
% % opts.StartPoint = [0.9, 1, 1, 1, 0.5]; 
% opts.StartPoint = [MaxEstimate pi/periodEstimate 5 0.011 0.53];
% % opts.Lower = [-Inf, -Inf, 0, 0, 0];
% % opts.Upper = [Inf, Inf, Inf, Inf, 1];
% opts.Lower = [0 0 0 0 -Inf 0];
% opts.Upper = [1 1000 Inf Inf Inf 1];
% 
% % Fit model to data
% [fitresult, gof] = fit(xData, yData, ft, opts);
% 
% end






