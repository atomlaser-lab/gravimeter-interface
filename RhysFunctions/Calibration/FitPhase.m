FigNum =20;
% 
y = r.data.R(:,1);
x = r.data.Param(1:numel(y));


x2 = linspace(0,180,1000);
[fitresult, ~] = createFit(x, y);

figure(FigNum);clf
scatter(x,y,'filled')
hold on
plot(x2,fitresult(x2))
xlabel('Phase')
ylabel('F = 2 pop')


Phi_midFringe = (fitresult.c + 3*pi/4)/pi*180/pi


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

