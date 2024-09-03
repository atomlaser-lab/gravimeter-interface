FigNum = 10;
% 
% y = r.data.R(:,1);
% x = r.data.Param(1:numel(y));

% y = r.data.C1.Nsum./(r.data.C1.Nsum + r.data.C2.Nsum);
% x = r.data.param(1:numel(y));

% y = r.data.raw.R(:,2,1,1);
% x = r.data.param1;
y = t0_7p5_T_10_w0_10.raw.R(:,end-2,1,1);
x = t0_7p5_T_10_w0_10.param1;

% y(5) = NaN;
% y(4) = NaN;
% y(end-6) = NaN;

x2 = linspace(0,260,1000);
[fitresult, ~] = createFit(x, y);
y_fit = fitresult(x2);

Phi_midFringe = (fitresult.c + 3*pi/4)/pi*180/pi*3/4

Min = x2(find(y_fit == min(y_fit),1));
Max = x2(find(y_fit == max(y_fit),1));

if Max > Min
    Phi_midFringe = (Max + Min)/2
else
    Phi_midFringe = (Min + Max)/2
end    


figure(FigNum);clf
scatter(x,y,'filled')
hold on
plot(x2,y_fit)
scatter(Phi_midFringe,fitresult(Phi_midFringe),'k','LineWidth',4)
xlabel('Phase')
ylabel('|p_0> pop')



annotation('textbox', [0.15, 0.1, 0.1, 0.1], 'String', sprintf('Mid-Fringe: %g deg, offet: %g, Contrast: %g, Phase Offset: %g', Phi_midFringe, fitresult.a, fitresult.b, fitresult.c ), ...
    'EdgeColor', 'none', 'BackgroundColor', 'white', 'FitBoxToText', 'on');
grid on
grid minor
if exist('TitleStuff','var') == 1
    Title = TitleStuff.Title;
        sgtitle({['{\bf\fontsize{14}' Title '}'],TitleStuff.SubTitle});
end

topaste(FigNum)

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

