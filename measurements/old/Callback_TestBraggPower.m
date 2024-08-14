function Callback_TestBraggPower(r)

if r.isInit()
    r.c.setup(500);
elseif r.isSet()
    
%     r.make();
%     r.upload;

    fprintf(1,'Run %d/%d\n',r.c.now,r.c.total);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    d = r.devices.d;
    d.getRAM;
    r.data.t(:,i1) = d.t;
    r.data.v(:,i1) = d.data(:,1);
    
    gauss = @(A,w,x0,x) A*exp(-(x-x0).^2/w^2);
    
%     nlf = nonlinfit(d.t,d.data(:,1),1e-3,d.t < 0.5e-4 || d.t > 1.5e-4);
%     nlf.setFitFunc(@(A,w,x0,y0,x) A*exp(-(x - x0).^2/w^2) + y0);
%     nlf.bounds2('A',[0,1,max(d.data(:,1)) - min(d.data(:,1))],'w',[10e-6,200e-6,30e-6],...
%         'x0',[0.5e-4,1.5e-4,1e-4],'y0',[-0.05,0.05,0]);
%     gauss = @(A,w,x0,x) A*exp(-(x-x0).^2/w^2);
%     nlf.setFitFunc(@(A1,A2,w,x0,f0,phi,y0,x) gauss(A1,w,x0,x) + gauss(A2,w,x0,x) + ...
%         2*sqrt(gauss(A1,w,x0,x).*gauss(A2,w,x0,x)).*sin(2*pi*f0*x + phi)+ y0);
%     nlf.bounds2('A1',[0,1,max(d.data(:,1)) - min(d.data(:,1))],'w',[10e-6,200e-6,30e-6],...
%         'x0',[0.5e-4,1.5e-4,1e-4],'y0',[-0.05,0.05,0],...
%         'A2',[0,1,max(d.data(:,1)) - min(d.data(:,1))],'phi',[-2*pi,2*pi,0],...
%         'f0',[2.4e5,2.75e5,2.57e5]);
    
    nlf = nonlinfit(d.t,d.data(:,1),1e-3,d.t < 0.5e-4 | d.t > 1.5e-4);
    nlf.setFitFunc(@(A1,w,x0,y0,x) gauss(A1,w,x0,x) +  y0);
    nlf.bounds2('A1',[0,1,max(d.data(:,1)) - min(d.data(:,1))],'w',[10e-6,200e-6,30e-6],...
        'x0',[0.5e-4,1.5e-4,1e-4],'y0',[-0.05,0.05,0]);
    nlf.fit,figure(1);clf;nlf.plot;
    r.data.c1(:,i1) = nlf.c(:,1);
    
    nlf.ex = d.t < 1.5e-4 | d.t > 3.25e-4;
    nlf.bounds2('A1',[0,1,max(d.data(:,1)) - min(d.data(:,1))],'w',[10e-6,200e-6,30e-6],...
        'x0',[2e-4,3e-4,2.5e-4],'y0',[-0.05,0.05,0]);
    nlf.fit,figure(2);clf;nlf.plot;
    r.data.c2(:,i1) = nlf.c(:,1);
    
    nlf.ex = d.t < 3.5e-4;

    nlf.setFitFunc(@(A1,A2,w,x0,f0,phi,y0,x) gauss(A1,w,x0,x) + gauss(A2,w,x0,x) + ...
        2*sqrt(gauss(A1,w,x0,x).*gauss(A2,w,x0,x)).*sin(2*pi*f0*x + phi)+ y0);
    nlf.bounds2('A1',[0,1,max(d.data(:,1)) - min(d.data(:,1))],'w',[10e-6,200e-6,30e-6],...
        'x0',[3.5e-4,5e-4,4e-4],'y0',[-0.05,0.05,0],...
        'A2',[0,1,max(d.data(:,1)) - min(d.data(:,1))],'phi',[-5*pi,5*pi,0],...
        'f0',[1.4e5,1.75e5,1.57e5]);

    nlf.fit,figure(3);clf;nlf.plot;
    
    r.data.c3(:,i1) = nlf.c(:,1);
    
    figure(10);clf;
%     subplot(1,2,1);
    plot(1:i1,r.data.c1(1,1:i1),'o-');
    hold on
    plot(1:i1,r.data.c2(1,1:i1),'sq-');
    plot(1:i1,r.data.c3([1,2],1:i1),'d-');
    ylim([0,Inf]);
    plot_format('Run','Amplitude [V]','',12);
%     subplot(1,2,2);
%     plot(1:i1,r.data.c(end-1,1:i1),'o-');
% %     ylim([0,Inf]);
%     plot_format('Run','Amplitude [V]','',12);
end