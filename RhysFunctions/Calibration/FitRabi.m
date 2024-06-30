


[x_sorted idx] = sort(PulseDuratrion.data.Param);
y_sorted = PulseDuratrion.data.R(idx,1).';
y_sorted(1:3) = 0.01;
y_sorted(y_sorted == 1) = nan;

[unique_x, ~, idx_uniq] = unique(x_sorted);
averaged_y = accumarray(idx_uniq, y_sorted, [], @mean);


[fitresult, gof] = createFit(unique_x, averaged_y);

x2 = linspace(unique_x(1),unique_x(end),numel(unique_x)*5);
y2 = fitresult(x2);


figure(10);clf
scatter(unique_x, averaged_y); % Example plot
hold on
plot(x2,y2)



function [fitresult, gof] = createFit(unique_x, averaged_y)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( unique_x, averaged_y );

% Set up fittype and options.
ft = fittype( 'a*sin(b*(x-g)+c)*exp(-d*(x-g))+f', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 -Inf];
opts.StartPoint = [0.859947605330542 0.823835382277067 0.405605888154915 0.0686301996375189 0.610717217585178 0.99098039626571];
opts.Upper = [1 2 Inf Inf 1 Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end
