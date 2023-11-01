function Callback_GenericOptimize(r)

if r.isInit()
    %Initialize run
    r.data.param = const.randomize(0.1:0.1:0.5);
    r.data.param2 = const.randomize(-25:2:-20);
    r.c.setup('var',r.data.param,r.data.param2);
    
elseif r.isSet()
    %Build/upload/run sequence
    r.make(r.devices.opt,'params',[r.data.param(r.c(1)),r.data.param2(r.c(2))]);
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, Param = %.3f, Param2 = %.3f\n',r.c.now,r.c.total,...
        r.data.param(r.c(1)),r.data.param2(r.c(2)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
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
    r.data.files{i1,i2} = img.raw.files;
    
    r.data.N(i1,i2) = img.get('N');
    r.data.F(i1,i2) = img.get('becFrac');
    r.data.T(i1,i2) = sqrt(prod(img.clouds.T));
    r.data.OD(i1,i2) = img.clouds.peakOD;
    
    figure(101);
    subplot(1,2,1);
    plot(r.data.param(1:i1),r.data.OD(1:i1,i2),'o');
    ylim([0,Inf]);
    xlim([min(r.data.param),max(r.data.param)]);
    grid on;
    
    if r.c.done(1)
        subplot(1,2,2);
        cla;
        str = {};
        for nn = 1:i2
            h = plot(r.data.param(1:i1),r.data.OD(1:i1,nn),'o');
            hold on;
            set(h,'MarkerFaceColor',h.Color);
            str{nn} = sprintf('%0.1f s',r.data.param2(nn));
        end
        legend(str);
        ylim([0,Inf]);
        xlim([min(r.data.param),max(r.data.param)]);
        grid on;
    end
end