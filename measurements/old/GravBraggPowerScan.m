function GravBraggPowerScan(r)

if r.isInit()
    r.data.param = const.randomize(0:0.001:0.015);
    r.c.setup('var',r.data.param);
elseif r.isSet()
    
    r.make(0,217.35e-3,1.35,r.data.param(r.c(1)),0,0,0);
    r.upload;
%     pause(5);
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,r.data.param(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5);
    c = Abs_Analysis_NClouds('last');
    
    r.data.N(i1,:) = [c.N];
    r.data.Nsum(i1,:) = [c.Nsum];
    r.data.R(i1,1) = r.data.N(i1,2)./sum(r.data.N(i1,:),2);
    r.data.R(i1,2) = r.data.Nsum(i1,2)./sum(r.data.Nsum(i1,:),2);
    
    figure(11);clf;
    subplot(1,2,1);
    plot(r.data.param(1:i1),r.data.N,'o');
    hold on
    plot(r.data.param(1:i1),r.data.Nsum,'sq');
    ylim([0,Inf]);
    subplot(1,2,2);
    plot(r.data.param(1:i1),r.data.R(:,1),'o');
    hold on;
    plot(r.data.param(1:i1),r.data.R(:,2),'sq');
    ylim([0,1]);
    grid on;
end