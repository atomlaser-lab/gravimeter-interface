function GravFringeMeasure(r)
%this function generates the spatial fringes and the spatial reference and
%measure the phase.


if r.isInit()
    
    alternating=zeros(20,1)+500e-6; %define a matrix 
    r.data.runcounter=linspace(1,length(alternating),length(alternating));
    r.data.param = alternating; %
    r.c.setup('var',r.data.param);
     %     clf(10);clf(11);clf(12);clf(15);clf(16);
clf(16);clf(15);clf(17)

elseif r.isSet()
    
    r.make(0,217e-3,1.078,0.176,r.data.param(r.c(1)));
    %for 217 ms drop time 
    r.upload;
%     pause(5);
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,1e3*r.data.param(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
%     pause(0.1);
    c = Abs_Analysis('last');
    f = c.fitdata;
    
    r.data.files{i1,1} = {c.raw.files(1).name,c.raw.files(2).name};
    r.data.c{i1,1} = c;
    r.data.imgs{i1,1} = c.imageCorr;
    r.data.y{i1,1} = c.fitdata.y;
    r.data.ydata{i1,1} = c.fitdata.ydata;
    
    
    r.data.ampref{i1,1}=c.fitdata.ydata(230-180:280-180); %y-axis data for the spatial reference! Careful to the ROI of the imaging system
    r.data.ampSF{i1,1}=c.fitdata.ydata(330-180:650-180); %y-axis data for the spatial fringes! Careful to the ROI of the imaging system
    
    y = f.ydata(abs(f.y-c.pos(2))<0.25e-3);
    r.data.contrast(i1,1) = (max(y) - min(y))./(max(y) + min(y) - 2*f.params.offset(2));
    
    figure(10);
    Ndim = [5,5];
    subplot(Ndim(1),Ndim(2),i1);
    c.plotAbsData([0,0.25],true);
    title(sprintf('%.3f ms',r.data.param(i1)*1e3));
    
%     figure(11);
%     subplot(Ndim(1),Ndim(2),i1);
%     plot(c.fitdata.y,c.fitdata.ydata,'.-');
%     title(sprintf('%.3f ms',r.data.param(i1)*1e3));
    
%     figure(12);clf;
%     plot(r.data.param(1:i1),r.data.contrast,'o-');
%     grid on;
    
    
%% Spatial fringes Fit

    %Spatial reference paremeters  
positref=230-180:1:280-180;    % "-180" value is due to the imaging ROI reseting the 0 reference.
height1=max(r.data.ampref{i1,1});
variance1=7;  %for the narrow spatial reference integrated OD its around 7
offset1=-010; %offset is around -17
peakpos1=83;  %peakposition depends on the ROI 

%Assymetric parameters
positSF=330-180:1:650-180;
ampsin0 = .76;
height0 = max(r.data.ampSF{i1,1})*ampsin0;
variance0 = 60;
offset0 = -18;
peakpos0 = 310;
freqsin0 = .23;
phasesin0 = -11;

    [V1SF, V2SF,SFpeakpos,enveloppepepos]= Integratedassymetricfit(positSF,r.data.ampSF{i1,1},...
        ampsin0,freqsin0,height0,offset0,peakpos0,phasesin0,variance0);
    
    [V1MZ, V2MZ,RefPeakpos]= Integratedspatialreffit(positref,r.data.ampref{i1,1},height1,peakpos1,variance1,offset1);
    
    r.data.FitREFrmse(i1,:) = V2MZ.rmse;
    r.data.SFrmse(i1,:) = V2SF.rmse ;
    r.data.enveloppepepos(i1,:)   =  enveloppepepos;
    r.data.RefPeakpos(i1,:) = RefPeakpos;
    r.data.SFpeakpos(i1,:) = SFpeakpos;
    r.data.delta_ref_SF(i1,:) = SFpeakpos-RefPeakpos;
    r.data.delta_env_ref(i1,:) = enveloppepepos-RefPeakpos;
    r.data.delta_sine_env(i1,:) = enveloppepepos-SFpeakpos;

%% Spatial fringes fit
pause(15e-3)
figure(15);
Ndim = [5,5];
subplot(Ndim(1),Ndim(2),i1);
% tiledlayout(5,2)
Integratedspatialreffit(positref,r.data.ampref{i1,1},height1,peakpos1,variance1,offset1);
hold off
pause(10e-3)
grid on;

figure(16);
Ndim = [5,5];
subplot(Ndim(1),Ndim(2),i1);
Integratedassymetricfit(positSF,r.data.ampSF{i1,1},ampsin0,freqsin0,height0,offset0,peakpos0,phasesin0,variance0);
hold off
pause(10e-3)
grid on;
%  
%  
% %% Fit details
% 
figure (17);
subplot(2,2,1)
plot(r.data.runcounter(1:i1),r.data.FitREFrmse,'*r',r.data.runcounter(1:i1),r.data.SFrmse,'ok');
% plot(r.data.param21(1:i1),r.data.FitREFrmse,'*r')
% hold on
% plot(r.data.SFrmse(1:i1),'ok')
%  h2= plot(r.data.param21(1:i1),SFrmse,'*r');
xlabel( 'Run ', 'Interpreter', 'none' );
ylabel( 'RMSE', 'Interpreter', 'none' );

xlim auto
ylim auto
grid on
% axis tight
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
hold on
legend('Reference','Spatial Fringes', 'Location', 'northeast');

hold off
title('Fits RMSE')

% subplot(1,1,1)
subplot(2,2,2)
% plot(r.data.param21,r.data.deltapeaks,'*r')
plot(r.data.runcounter(1:1:i1),r.data.delta_ref_SF,'*r')
% hold on
xlabel( 'Run ', 'Interpreter', 'none' );
ylabel( 'Difference in pixels', 'Interpreter', 'none' );
xlim auto
ylim auto
grid on
%  axis tight
 curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
title('\Delta difference in position between spatial reference and sine modulation peak')
hold off

subplot(2,2,3)
% hold off
plot(r.data.runcounter(1:1:i1),r.data.delta_env_ref,'*r')
% hold off
xlabel( 'Run ', 'Interpreter', 'none' );
ylabel( 'Difference in pixels', 'Interpreter', 'none' );
xlim auto
ylim auto
grid on
%  axis tight
 curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
title('\Delta difference in position between spatial reference and SF gaussian')

subplot(2,2,4)
plot(r.data.runcounter(1:1:i1),r.data.delta_sine_env,'*r')
% hold on
xlabel( 'Run ', 'Interpreter', 'none' );
ylabel( 'Difference in pixels', 'Interpreter', 'none' );
xlim auto
ylim auto
grid on
%  axis tight
 curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
title('\Delta difference in position between SF sine peak and SF gaussian')

   
end