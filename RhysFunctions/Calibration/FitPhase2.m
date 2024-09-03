FigNum = 10;

ClearFigOnOff = 0;
% data = t0_0_T_10_w0_10;
% data = t0_7p5_T_10_w0_10;
data = t0_15_T_10_w0_10;


XX = data.param2;
y = data.raw.R(:,:,1,1);
x = data.param1;
x2 = linspace(0,200,100);

N_Total = data.raw.N(:,1,1,2) + data.raw.N(:,1,1,1);
Nsum_Total = data.raw.Nsum(:,1,1,2) + data.raw.Nsum(:,1,1,1);

for ii = 1:size(y,2)
% % % Remove outliers
%     Remove N = 0



    warning('off','all')
    [fitresult, ~] = createFit(x, y(:,ii,1,1));
    warning('on','all')

    ci = confint(fitresult);
    a(ii) = fitresult.a;
    b(ii) = fitresult.b;
    c(ii) = fitresult.c;

    a_uncertainty(ii) = (ci(2,1) - ci(1,1)) / 2;
    b_uncertainty(ii) = (ci(2,2) - ci(1,2)) / 2;
    c_uncertainty(ii) = (ci(2,3) - ci(1,3)) / 2; 

%     y_fit(:,ii) = fitresult(x2);
end

figure(FigNum);
if ClearFigOnOff == 1
    clf
end
%%
subplot(3,1,1)
errorbar(XX,a,a_uncertainty)
hold on
scatter(XX,a)
grid on
grid minor
ylabel('Population Offset','Interpreter','latex','FontSize',14)
ylim([0 0.5])

subplot(3,1,2)
errorbar(XX,b,b_uncertainty)
hold on
scatter(XX,b)
grid on
grid minor
ylabel('Contrast','Interpreter','latex','FontSize',14)
ylim([0,1])

subplot(3,1,3)
errorbar(XX,c,c_uncertainty)
hold on
scatter(XX,c)
grid on
grid minor
ylabel('Phase Offset','Interpreter','latex','FontSize',14)

xlabel('r_final /w_0','Interpreter','latex','FontSize',14)


clear a b c a_uncertainty b_uncertainty c_uncertainty




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

