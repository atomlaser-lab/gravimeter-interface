function Callback_MeasureMWpower_update(r)

if r.isInit()
    
    r.data.param = const.randomize(100:100:700);
r.data.ParamName = 'Pulse duration';
r.data.ParamUnits = 'us';
    r.c.setup('var', r.data.param);
elseif r.isSet()
    
    r.make(r.devices.opt,'params',r.data.param(r.c(1))).upload;
   
fprintf(1,'Run %d/%d, %s = %.3f%s \n ',r.c.now,r.c.total,...
        r.data.ParamName,r.data.param(r.c(1)),r.data.ParamUnits);


elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.25 + 0.1*rand);

    img = Abs_Analysis_GUI('last',1);

    
    r.data.files{i1,1} = img.raw.files;
    
%     Get processed data
    
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));

    figure(98);clf;
    subplot(1,2,1)
    plot(r.data.param(1:i1),r.data.R(1:i1,:),'o');
    enhformat('Pulse duration [us]','Population');
% %     h = legend('m = -1','m = 0');
%     set(h,'Location','West');
    ylim([0,1]);
    title(' Microwave frequency using fit over OD')
    grid on
    hold on;
    
%     subplot(1,3,2)
%     plot(r.data.df(1:i1),r.data.R(1:i1,:),'sq');
%     hold off;
%     plot_format('Freq [MHz]','Population','',12);
%     h = legend('m = -1','m = 0');
%     set(h,'Location','West');
%     title(' Raman frequency using ROI')
%     grid on;
%     if r.c.done
%     tNow = datestr(now);
%         caption = sprintf('Determination of Raman frequency %s', tNow);
%         sgtitle(caption)
%     end
    
    subplot(1,2,2)
    plot(r.data.param(1:i1),r.data.N(1:i1,:),'o');
    enhformat('Freq [kHz]',' Atom Number');
    h = legend('m = -1','m = 0');
    set(h,'Location','West');
    grid on
    hold on;
    
%     if r.c.done(1)
%         nlf = nonlinfit(r.data.df*1e-3,r.data.R(:,1),1e-2);
%         nlf.setFitFunc(@(A,R,x0,x) A*(1 - 4*R.^2./(4*R^2+(x-x0).^2).*sin(2*pi*sqrt(4*R^2+(x-x0).^2).*372/2).^2));
%         [~,idx] = min(nlf.y);
%         nlf.bounds2('A',[0.9,1,0.95],'R',[0,10,0.8]*1e-3,'x0',[min(nlf.x),max(nlf.x),nlf.x(idx)]);
%         nlf.fit;
%         fprintf(1,'Rabi frequency: %.3f kHz, Center = %.3f kHz\n',nlf.c(2,1)*1e3,nlf.c(3,1)*1e3);
%         subplot(1,2,1);
%         hold on
%         xplot = linspace(min(nlf.x),max(nlf.x),1e2);
%         plot(xplot*1e3,nlf.f(xplot),'-');
%         hold off;
% % 
% %         nlf = nonlinfit(r.data.freq2 - const.f_Rb_groundHFS/1e6,r.data.N/max(r.data.N));
% %         nlf.setFitFunc(@(A,R,x0,x) A*4*R.^2./(4*R^2+(x-x0).^2).*sin(2*pi*sqrt(4*R^2+(x-x0).^2).*215/2).^2);
% %         [~,idx] = max(nlf.y);
% %         nlf.bounds2('A',[0.9,1,0.95],'R',[0,10,0.8]*1e-3,'x0',[-0.003,0.003,nlf.x(idx)]);
% %         nlf.fit;
% %         fprintf(1,'Rabi frequency: %.3f kHz, Center = %.3f kHz\n',nlf.c(2,1)*1e3,nlf.c(3,1)*1e3);
% %         figure(1);clf;
% %         nlf.plot;
%         
%     end
    
%     Ndim = ceil(sqrt(r.c.total));
%     figure(29);
%     subplot(Ndim,Ndim,i1);

%     title(sprintf('Chirp srate = %.3f MHz',r.data.chirp(i1)/1e6));
end