function Callback_MeasureImagingFrequency(r)

if r.isInit()
    %Initialize run
%     r.data.detuning = const.randomize(unique([-16:2:16,-6:1:6]));
        r.data.detuning = const.randomize(-3:.2:3);

    r.c.setup('var',r.data.detuning);
    
%     r.makerCallback = @makeSequence;
elseif r.isSet()
    %Build/upload/run sequence
    r.make(r.devices.opt,'detuning',r.data.detuning(r.c(1)));
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, Detuning = %.1f MHz\n',r.c.now,r.c.total,...
        r.data.detuning(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_GUI('last');
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
    r.data.files{i1,1} = img.raw.files;
    
    r.data.N(i1,1) = img.get('N');
    r.data.T(i1,1) = sqrt(prod(img.clouds.T));
    r.data.OD(i1,1) = img.clouds.peakOD;
    
    figure(101);clf;
    subplot(1,2,1)
    plot(r.data.detuning(1:i1),r.data.N(1:i1)/1e6,'o');
    ylim([0,Inf]);
    xlim([min(r.data.detuning),max(r.data.detuning)]);
    grid on;
    enhformat('Detuning [MHz]','Number [x$10^6$]')
    
 subplot(1,2,2)
    plot(r.data.detuning(1:i1),r.data.OD(1:i1),'o');
    ylim([0,Inf]);
    xlim([min(r.data.detuning),max(r.data.detuning)]);
    grid on;
    enhformat('Detuning [MHz]','OD')


%     if r.c.done(1) || i1 > 5
%         subplot(1,2,1)
%         nlf = nonlinfit(r.data.detuning(1:i1),r.data.N(1:i1)/1e6,0.05,r.data.N == 0 & abs(r.data.detuning(1:i1)) < 5);
%         nlf.setFitFunc(@(A,w,x0,x) A./(1 + 4*(x-x0).^2/w^2));
%         nlf.bounds2('A',[0,2*max(nlf.y),max(nlf.y)],'w',[0,20,6],'x0',[min(nlf.x),max(nlf.x),0]);
%         nlf.fit
%         hold on;
%         xx = linspace(min(nlf.x),max(nlf.x),100);
%         plot(xx,nlf.f(xx),'--','linewidth',2);
%     end

  if r.c.now == r.c.total
        subplot(1,2,1)
        nlf = nonlinfit(r.data.detuning(1:i1),r.data.N(1:i1)/1e6,0.05,r.data.N == 0 & abs(r.data.detuning(1:i1)) < 5);
        nlf.setFitFunc(@(A,w,x0,x) A./(1 + 4*(x-x0).^2/w^2));
        nlf.bounds2('A',[0,2*max(nlf.y),max(nlf.y)],'w',[0,20,6],'x0',[min(nlf.x),max(nlf.x),0]);
        nlf.fit
        hold on;
        xx = linspace(min(nlf.x),max(nlf.x),100);
        plot(xx,nlf.f(xx),'--','linewidth',2);
    end

end