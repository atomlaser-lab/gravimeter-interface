function Callback_MeasureTrapFrequency(r)

if r.isInit()
    %Initialize run
%     r.data.param = [2.5e-3:2.5e-3:50e-3];
    r.data.time = 1e-3:2e-3:40e-3;
    r.c.setup('var',r.data.time);
    
    r.makerCallback = @makeSequence;
elseif r.isSet()
    %Build/upload/run sequence
%     r.make(0,216.5e-3,1.48,r.data.param(r.c(1)));
    r.make(r.devices.opt.set('params',r.data.time(r.c(1))));
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, Delay = %.3f ms\n',r.c.now,r.c.total,...
        r.data.time(r.c(1))*1e3);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis('last');
    if ~img.raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif r.c.now > 1 && strcmpi(r.data.files{r.c.now - 1}.name,img.raw.files.name)
        r.c.decrement;
        return
    end
    
    r.data.w(i1,:) = img.clouds(1).becWidth;
    r.data.wg(i1,:) = img.clouds(1).gaussWidth;
    r.data.pos(i1,:) = img.clouds(1).pos;
    r.data.files{i1,1} = img.raw.files;
    
    figure(11);clf;
    subplot(1,2,1);
    plot(r.data.time(1:i1)*1e3,1e6*r.data.w(1:i1,:),'o-');
    grid on;
    plot_format('Delay [ms]','Width [um]','',12);
    legend('x','y');
    ylim([0,1000]);
    subplot(1,2,2);
    plot(r.data.time(1:i1)*1e3,1e6*(r.data.pos(1:i1,:) - r.data.pos(1,:)),'o-');
    grid on;
    plot_format('Delay [ms]','Position [um]','',12);
    legend('x','y');
%     ylim([0,1000]);

end