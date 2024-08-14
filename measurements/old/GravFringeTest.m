function GravFringeTest(r)

if r.isInit()
    
    alternating=zeros(10,1);
    alternating(1:2:end)=500e-6; %put 500e-6 asymetry every second run in the matrix
    runcounter=linspace(1,length(alternating),length(alternating));
    
    r.data.param = alternating;
    r.data.param21=runcounter;
    r.c.setup('var',r.data.param);
    
    clf(10);clf(11);clf(12);clf(15);clf(16);
elseif r.isSet()
    
    r.make(0,216.5e-3,1.075,0.176,r.data.param(r.c(1)));
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
    
    r.data.amp{i1,1}=c.imageCorr(150:500,715);
    y = f.ydata(abs(f.y-c.pos(2))<0.25e-3);
    r.data.contrast(i1,1) = (max(y) - min(y))./(max(y) + min(y) - 2*f.params.offset(2));
    
    figure(10);
    Ndim = [5,3];
    subplot(Ndim(1),Ndim(2),i1);
    c.plotAbsData([0,0.25],true);
    title(sprintf('%.3f ms',r.data.param(i1)*1e3));
    
    figure(11);
    subplot(Ndim(1),Ndim(2),i1);
    plot(c.fitdata.y,c.fitdata.ydata,'.-');
    title(sprintf('%.3f ms',r.data.param(i1)*1e3));
%     
    figure(12);clf;
    plot(r.data.param(1:i1),r.data.contrast,'o-');
    grid on;
    
    
%%Spatial fringes Fit
    %Synetric parameters   
posit=0:1:350;    
xvec=posit;
height1=max(r.data.amp{i1,1});
variance1=150;
offset1=-0.05;
peakpos1=170;
%Assymetric parameters
xvec = posit;
height0 = max(r.data.amp{i1,1});
variance0 = 72.21;
offset0 = -0.059;
peakpos0 = 170.1;
freqsin0 = 0.2363;
ampsin0 = -0.534;
phasesin0 = 6.355;
if mod(r.c(1),2)==0
    [V1MZ, V2MZ,PeakposMZ]=symetricfit(posit,r.data.amp{i1,1},height1,peakpos1,variance1,offset1);
    r.data.Fitrmse=V2MZ.rmse;
else
    [V1SF, V2SF,SFpeakpos,enveloppepepos]= assymetricfit(posit,r.data.amp{i1,1},ampsin0,freqsin0,height0,offset0,peakpos0,phasesin0,variance0);
    r.data.SFrmse=V2SF.rmse;
end


figure(15);
Ndim = [5,2];
subplot(Ndim(1),Ndim(2),i1);
 if mod(r.c(1),2)==0
symetricfit(posit,r.data.amp{i1,1},height1,peakpos1,variance1,offset1);
pause(10e-3)
grid on;
 else
assymetricfit(posit,r.data.amp{i1,1},ampsin0,freqsin0,height0,offset0,peakpos0,phasesin0,variance0);
pause(10e-3)
grid on
 end
 
 
%% Fit details
figure (16);
% subplot(2,2,[1,2])
subplot(1,1,1)
if mod(r.c(1),2)==0
   h1= plot(r.data.param21(2:2:i1),r.data.Fitrmse,'ok');
else
    hold off
  h2= plot(r.data.param21(1:2:i1),r.data.SFrmse,'*r');
end
xlabel( 'Run ', 'Interpreter', 'none' );
ylabel( 'RMSE', 'Interpreter', 'none' );
xlim auto
ylim auto
grid on
% axis tight
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
hold on
%just to define a f*****g legend
LH(1) = plot(nan, nan, '*r');
L{1} = 'MZ';
LH(2) = plot(nan, nan, 'ok');
L{2} = 'SF';
legend(LH,L, 'Location', 'northeast');
title('Fits RMSE')

% subplot(2,2,3)
% 
% plot(plot(r.data.param21(1:2:i1),)
% title('\Delta_SF difference in position between gaussian and sine')

% subplot(2,2,4)
% plot(r.c(1),enveloppepeakpos-PeakposMZ)
% title('\Delta difference in position between MZgaussian and SFgaussian')

   
end