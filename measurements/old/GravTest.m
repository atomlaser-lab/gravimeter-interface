function GravTest(r)

if r.isInit()
    r.data.param = const.randomize(0.1:0.05:0.5);
    r.c.setup('var',r.data.param);
elseif r.isSet()
    
    r.make(0,216.5e-3,1.05,r.data.param(r.c(1)));
    r.upload;
%     pause(5);
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,r.data.param(r.c(1)));
    
elseif r.isAnalyze()
    nn = r.c(1);
    pause(0.5);
    c = Abs_Analysis('last');
%     r.data.c(nn,1) = c;
%     r.data.x(:,nn) = c.fitdata.x;
%     r.data.xdata(:,nn) = c.fitdata.xdata;
    
    r.data.N(nn,:) = c.N;
    r.data.Nsum(nn,:) = c.Nsum;
%     r.data.R(nn,1) = r.data.N(nn,2)./sum(r.data.N(nn,:),2);
%     r.data.Nbec(nn,1) = c.N.*c.becFrac;
    r.data.T(nn,:) = c.T;
    r.data.OD(nn,1) = c.peakOD;
%     r.data.pos(nn,:) = c.pos;
    r.data.PSD(nn,1) = c.PSD;
    
%     figure(10);%clf;
%     Nsub = ceil(sqrt(r.c.total));
% %     for mm = 1:nn
%         subplot(Nsub,Nsub,nn);
% %         plot(r.data.y(:,nn),r.data.ydata(:,nn),'.-');
%         c.plotAbsData([0,5],true);
%         title(sprintf('T = %.1f',r.data.param(r.c.now)*1e6));
%     end
    
    figure(11);clf;
    subplot(1,3,1);
    plot(r.data.param(1:nn),r.data.N(1:nn),'o');
    ylim([0,Inf]);
    subplot(1,3,2);
    plot(r.data.param(1:nn),r.data.T(1:nn,:)*1e6,'o');
    ylim([0,Inf]);
    subplot(1,3,3);
    plot(r.data.param(1:nn),r.data.PSD(1:nn),'o');
    ylim([0,Inf]);
end