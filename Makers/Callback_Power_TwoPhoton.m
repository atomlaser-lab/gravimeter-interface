function Callback_Power_TwoPhoton(r)

if r.isInit()

    r.data.Power = [0.1:0.1:1]; %for dipoles: [1.5:0.5:8] %for RF: [0.5:0.25:2]
    r.data.delta = const.randomize(linspace(19.95,20.1,20));% phase AOM %19
    r.data.time = const.randomize(linspace(0,96,25))*1e-6;
    r.data.param = [1:length(r.data.delta) + length(r.data.time)];
    r.data.limit = 1;
    r.c.setup('var',r.data.param,r.data.Power);

elseif r.isSet()

    if r.c(1) < length(r.data.delta) + 1
        r.make(r.devices.opt,'param1',r.data.delta(r.c(1)),'param2',40e-6,'param3',r.data.Power(r.c(2))).upload;
        fprintf(1,'Run %d/%d, delta = %.3f MHz, Pulse Time = %.3f microsec, P = %.1f\n',r.c.now,r.c.total,r.data.delta(r.c(1)),40,r.data.Power(r.c(2)));
    else
        r.make(r.devices.opt,'param1',r.data.two_photon_opt(r.c(2)),'param2',r.data.time(r.c(1) - r.data.counter),'param3',r.data.Power(r.c(2))).upload;
        fprintf(1,'Run %d/%d, delta = %.3f MHz, Pulse Time = %.0f microsec, P = %.1f\n',r.c.now,r.c.total,r.data.two_photon_opt(r.c(2)),1e6*r.data.time(r.c(1) - r.data.counter),r.data.Power(r.c(2)));
    end

elseif r.isAnalyze()

    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_DualState('last',1);

    r.data.files{i1,i2,1} = img(1).raw.files;
    r.data.N(i1,i2,:) = img.get('N');
    r.data.Nsum(i1,i2,:) = img.get('Nsum');
    r.data.peakOD(i1,i2,:) = img.get('peakOD');
    r.data.R(i1,i2,:) = r.data.N(i1,i2,:)./sum(r.data.N(i1,i2,:));
    r.data.Rsum(i1,i2,:) = r.data.Nsum(i1,i2,:)./sum(r.data.Nsum(i1,i2,:));

    if (~img(1).raw.status.ok() || r.data.N(i1,i2,1) == 0 || r.data.N(i1,i2,1) > 1e7 || r.data.Rsum(i1,i2,1) > 0.98 || r.data.Rsum(i1,i2,2) > 0.98)
        %         if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0)
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        if (~img(1).raw.status.ok())
            warning('Imaging failed!')
            r.c.decrement;

        elseif (r.data.N(i1,1) == 0 || r.data.N(i1,i2,1) > 1e7 || r.data.Rsum(i1,i2,1) > 0.98 || r.data.Rsum(i1,i2,2) > 0.98)
            %             elseif (r.data.N(i1,1) == 0)
            r.data.limit=1+r.data.limit;
            warning('Number of atoms too low!')
            r.c.decrement;
            if r.data.limit == 50
                r.data.limit = 0;
                r.stop
                error('Run failed! lock your lasers.')
            end
        end
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end

    if r.c(1) < length(r.data.delta) + 1

        figure(42)
        clf
        hold on
        scatter(r.data.delta(1:i1),r.data.Rsum([1:i1],i2,2),'filled')
        xlim([min(r.data.delta) max(r.data.delta)])
        ylim([0 1])
        xlabel('Two Photon Detuning (MHz)')
        ylabel('N_2')
        title(sprintf('LG Power = %.1f',r.data.Power(i2)))

        r.data.counter = i1;
        r.data.two_photon_opt(i2) = r.data.delta(find(r.data.Rsum([1:i1],i2,2) == max(r.data.Rsum([1:i1],i2,2))));

    else

        figure(43)
        clf
        hold on
        scatter(r.data.time(1:r.c(1) - r.data.counter),r.data.Rsum([1:r.c(1) - r.data.counter],i2,2),'filled')
        xlim([min(r.data.time) max(r.data.time)])
        ylim([0 1])
        xlabel('Pulse Time (microsec)')
        ylabel('N_2')
        title(sprintf('LG Power = %.1f',r.data.Power(i2)))

%         r.data.pulse_time_opt(i1,i2) = r.data.time(find(r.data.Rsum([1:i1-r.data.counter],i2,2) == max(r.data.Rsum([1:i1],i2,2))));

    end

    if r.c(1) == length(r.data.param)

        figure(42)
        saveas(gcf,sprintf('D:\\data\\VMG_autosave\\P =  %.2f, Two Photon.png',r.data.Power(i2)));
        figure(43)
        saveas(gcf,sprintf('D:\\data\\VMG_autosave\\P =  %.2f, Pulse Time.png',r.data.Power(i2)));

    end

    if r.c.done(1)

        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\data.mat','data');
        saveas(gcf,'D:\data\VMG_autosave\data.fig');

    end

end