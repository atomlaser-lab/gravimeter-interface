
x = TwoPhotonScan_2.data.Param.';
y_av = (TwoPhotonScan_1.data.R(:,1) + TwoPhotonScan_2.data.R(:,1))/2;

[fitresult, gof] = createFit(x, y_av);

x2 = linspace(x(1),x(end),numel(x)*10);
y2 = fitresult(x2);

figure(10);clf
scatter(x,y_av)
hold on
plot(x2,y2)
xlabel('Detuning (kHz)','Interpreter','latex','FontSize',16)
ylabel('F = 2 Population','Interpreter','latex','FontSize',16)





function [fitresult, gof] = createFit(x, y)

% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a*sinc(b*(x+d))^2+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [1 0 0 -10];
opts.StartPoint = [1 12 0 0.877484639967382];
opts.Upper = [1 100 0 10];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


end

