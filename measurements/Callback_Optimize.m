function Callback_Optimize(r)

if r.isInit()
    r.data.tof = 30e-3;
    r.data.param1 = const.randomize(0:5:40);
    r.data.param2 = const.randomize(0:0.5:2);
    r.c.setup('var',r.data.param1,r.data.param2);
elseif r.isSet()
    r.make(r.data.tof,r.data.param1(r.c(1)),r.data.param2(r.c(2))).upload;
    fprintf(1,'Run %d/%d, Param1 = %.3f, Param2 = %.3f\n',r.c.now,r.c.total,r.data.param1(r.c(1)),r.data.param2(r.c(2)));
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.25);
    img = Abs_Analysis('last',1);
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end

    r.data.files{i1,i2} = img.raw.files;
    r.data.N(i1,i2) = img.get('N');
    r.data.T(i1,i2) = sqrt(prod(squeeze(img.get('T'))));
    r.data.PSD(i1,i2) = img.get('PSD');
    r.data.OD(i1,i2) = img.get('peakOD');

    figure(123);
    if r.c.now() == 1
        clf;
    end
    subplot(2,2,[1,3]);
    errorbar(r.data.param1(1:i1),r.data.N(1:i1,i2),0.05*r.data.N(1:i1,i2),'o');
    plot_format('Param 1','Number','',12);
    ylim([0,Inf]);
    grid on;

    if r.c.done(1)
        subplot(2,2,2);
        cla;
        s = {};
        for nn = 1:i2
            errorbar(r.data.param1(1:i1),r.data.OD(1:i1,nn),0.05*r.data.OD(1:i1,nn),'o');
            s{nn} = sprintf('Param 2 = %.2f',r.data.param2(nn));
            hold on
        end
        plot_format('Param 1','OD','',12);
        ylim([0,Inf]);
        grid on;
        legend(s);

        subplot(2,2,4);
        cla;
        s = {};
        for nn = 1:i2
            errorbar(r.data.param1(1:i1),r.data.N(1:i1,nn),0.05*r.data.N(1:i1,nn),'o');
            s{nn} = sprintf('Param 2 = %.2f',r.data.param2(nn));
            hold on
        end
        plot_format('Param 1','N','',12);
        ylim([0,Inf]);
        grid on;
        legend(s);
    end

%     figure(2);clf;
%     plot(r.data.param1(1:i1),r.data.PSD(1:i1,:),'o');
%     ylim([0,Inf]);
%     grid on;
end


end