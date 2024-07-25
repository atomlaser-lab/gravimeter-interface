% % % Inputs
MovingMean = 20;


% % % Contrast Data
x = T_100us.data.Param;
y = T_100us.data.R(:,1);

[fitresult, ~] = ContrastFit(x,y);

Contrast = fitresult.b;
Offset = fitresult.a;

x2 = linspace(0,180,1000);

figure(11);clf
subplot(2,1,1)
scatter(x,y)
hold on
plot(x2,fitresult(x2))


%%
% % % Al. Dev Data
N1 = T_100_stability.data.R(:,1);
N2 = T_100_stability.data.R(:,1);

y2 =  T_100_stability.data.N(:,1);
y3 =  T_100_stability.data.N(:,2);
x = T_100_stability.data.Param(1:numel(N1));

% % % remove bad data points
figure(10);clf
subplot(2,1,1);hold on
scatter(x(2:end),abs(diff(y3)),'b')


Mean_dif = mean(abs(diff(y3)));
std_dif = std(abs(diff(y3)));

BadRunPos = find(abs(diff(y3)) > Mean_dif + 3*std_dif );

y3(BadRunPos) = NaN;
y2(BadRunPos) = NaN;
N1(BadRunPos) = NaN;
N2(BadRunPos) = NaN;

scatter(x(2:end),abs(diff(y3)),'r')
legend('All Runs','Bad Runs Removed')


% % % Plot Data
subplot(2,1,2); hold on
scatter(x,N1,'r')
plot(x,movmean(N1,MovingMean),'b','LineWidth',2)
legend('Raw','Moving Mean')
xlabel('Run')
ylabel('F = 2 Population')

% % % Convert Pop to Phase
phase = acos((2*N1 - Offset)/Contrast);

figure(11);
subplot(2,1,2)
scatter(x,phase,'r')
plot(x,movmean(phase,MovingMean),'b','LineWidth',2)
legend('Raw','Moving Mean')
xlabel('Run')
ylabel('Phase')






%% Contrast Fit
function [fitresult, gof] = ContrastFit(x, y)
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a+b*sin(x*pi/180+c)^2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -3.14];
opts.StartPoint = [0.408719846112552 0.594896074008614 0.262211747780845];
opts.Upper = [1 1 3.14];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end


