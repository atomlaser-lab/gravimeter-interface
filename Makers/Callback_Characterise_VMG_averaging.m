function Callback_Characterise_VMG_averaging(r)

if r.isInit()

    r.data.no_averages = 2;
    r.data.delta_param = 2;
    r.data.delta = repmat(linspace(19.9,20.15,r.data.delta_param),1,r.data.no_averages);
    r.data.delta_key = randperm(r.data.delta_param*r.data.no_averages); %getting random key
    r.data.delta = r.data.delta(r.data.delta_key); %randomising using key
    r.data.time_param = 26;
    r.data.time = repmat(linspace(32,32,r.data.time_param)*1e-6,1,r.data.no_averages);
    r.data.time_key = randperm(r.data.time_param*r.data.no_averages); %getting random key
    r.data.time = r.data.time(r.data.time_key); %randomising using key
    r.data.phase_param = 5;
    r.data.phase = repmat(linspace(0,180,r.data.phase_param),1,r.data.no_averages);
    r.data.phase_key = randperm(r.data.phase_param*r.data.no_averages); %getting random key
    r.data.phase = r.data.phase(r.data.phase_key); %randomising using key
    r.data.limit = 1; %error limit
    r.data.param = [1:(r.data.delta_param+r.data.time_param+r.data.phase_param)*r.data.no_averages];
    r.c.setup('var',r.data.param);

%     r.data.no_averages = 10;
%     r.data.delta_param = 26;
%     r.data.delta = repmat(linspace(19.9,20.15,r.data.delta_param),1,r.data.no_averages);
%     r.data.delta_key = randperm(r.data.delta_param*r.data.no_averages); %getting random key
%     r.data.time_param = 26;
%     r.data.time = repmat(linspace(0,100,r.data.time_param)*1e-6,1,r.data.no_averages);
%     r.data.time_key = randperm(r.data.time_param*r.data.no_averages); %getting random key
%     r.data.phase_param = 5;
%     r.data.phase = repmat(linspace(0,180,r.data.phase_param),1,r.data.no_averages);
%     r.data.phase_key = randperm(r.data.phase_param*r.data.no_averages); %getting random key
%     r.data.limit = 1; %error limit
%     r.data.phase = r.data.phase(r.data.phase_key);
%     r.data.param = [1:(r.data.delta_param+r.data.time_param+r.data.phase_param)*r.data.no_averages];
%     r.c.setup('var',r.data.param);

elseif r.isSet()

    if r.c(1) < r.data.no_averages*r.data.delta_param+1
        r.make(r.devices.opt,'param1',r.data.delta(r.c(1)),'param2',36e-6,'param3',0,'params',1).upload;
        fprintf(1,'Run %d/%d, delta = %.3f MHz, time = %.1f microsec, phase = %.3f degrees\n',r.c.now,r.c.total,r.data.delta(r.c(1)),36,0);
    elseif ((r.c(1) > r.data.no_averages*r.data.delta_param) & (r.c(1) < r.data.no_averages*(r.data.delta_param + r.data.time_param) + 1))
        r.make(r.devices.opt,'param1',r.data.two_photon_opt,'param2',r.data.time(r.c(1) - r.data.counter),'param3',0,'params',1).upload;
        fprintf(1,'Run %d/%d, delta = %.3f MHz, time = %.1f microsec, phase = %.3f degrees\n',r.c.now,r.c.total,r.data.two_photon_opt,r.data.time(r.c(1) - r.data.counter)*1e6,0);
    elseif ((r.c(1) > r.data.no_averages*(r.data.delta_param + r.data.time_param)) & (r.c(1) < r.data.no_averages*(r.data.delta_param + r.data.time_param + r.data.phase_param) + 1))
        r.make(r.devices.opt,'param1',r.data.two_photon_opt,'param2',r.data.time_opt,'param3',r.data.phase(r.c(1) - r.data.counter),'params',0).upload;
        fprintf(1,'Run %d/%d, delta = %.3f MHz, time = %.1f microsec, phase = %.3f degrees\n',r.c.now,r.c.total,r.data.two_photon_opt,1e6*r.data.time_opt,r.data.phase(r.c(1) - r.data.counter));
    end

elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_DualState('last',1);

    r.data.files{i1,1} = img(1).raw.files;
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.peakOD(i1,:) = img.get('peakOD');
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));

%     if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0 || r.data.N(i1,1) > 2e7)
%         %         if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0)
%         %
%         % Checks for an error in loading the files (caused by a missed
%         % image) and reruns the last sequence
%         %
%         if (~img(1).raw.status.ok())
%             warning('Imaging failed!')
%             r.c.decrement;
% 
%         elseif r.data.N(i1,1) == 0 || r.data.N(i1,1) > 2e7
%             r.data.limit=1+r.data.limit;
%             warning('Number of atoms too low!')
%             r.c.decrement;
%             if r.data.limit == 20
%                 r.data.limit = 0;
%                 r.stop
%                 error('Run failed! lock your lasers.')
%             end
%         end
%         return;
%     elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
%         r.c.decrement;
%         return;
%     end

    if i1 == 1
        r.data.N1_delta = NaN*ones(r.data.delta_param,r.data.no_averages);
        r.data.N2_delta = NaN*ones(r.data.delta_param,r.data.no_averages);
        r.data.N1_time = NaN*ones(r.data.time_param,r.data.no_averages);
        r.data.N2_time = NaN*ones(r.data.time_param,r.data.no_averages);
        r.data.N1_phase = NaN*ones(r.data.phase_param,r.data.no_averages);
        r.data.N2_phase = NaN*ones(r.data.phase_param,r.data.no_averages);
    end

    if r.c(1) < r.data.no_averages*r.data.delta_param+1
        % plotting as you go
        A1=r.data.delta_key';
        B1=r.data.Rsum(:,2);
        [A1,index(1:i1)] = sortrows(A1(1:i1));
        B1 = B1(index(i1),:);
        A2=r.data.delta_key';
        B2=r.data.Rsum(:,1);
        [A2,index(1:i1)] = sortrows(A2(1:i1));
        B2 = B2(index(i1),:);
        x = linspace(min(r.data.delta),max(r.data.delta),r.data.delta_param);

        figure(41)
        clf
        hold on
        r.data.N1_delta(r.data.delta_key(i1)) = r.data.Rsum(i1,1);
        r.data.N2_delta(r.data.delta_key(i1)) = r.data.Rsum(i1,2);
        C2_mean = mean(r.data.N2_delta,2,"omitnan");
        C2_std = std(r.data.N2_delta,0,2,"omitnan");
        errorbar(x,C2_mean,C2_std,"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
        xlim([min(linspace(19.9,20.15,r.data.delta_param)) max(linspace(19.9,20.15,r.data.delta_param))])
        ylim([0 1])
        xlabel('Two Photon Detuning (MHz)')
        ylabel('N_2')

    elseif ((r.c(1) > r.data.no_averages*r.data.delta_param) & (r.c(1) < r.data.no_averages*(r.data.delta_param + r.data.time_param) + 1))

        counter = i1 - r.data.no_averages*r.data.time_param;

        % plotting as you go
        A1=r.data.time_key';
        B1=r.data.Rsum([r.data.no_averages*r.data.time_param:i1],2);
        [A1,index(1:counter)] = sortrows(A1(1:counter));
        B1 = B1(index(counter),:);
        A2=r.data.time_key';
        B2=r.data.Rsum([r.data.no_averages*r.data.time_param:i1],1);
        [A2,index(1:counter)] = sortrows(A2(1:counter));
        B2 = B2(index(counter),:);
        x = linspace(min(r.data.time),max(r.data.time),r.data.time_param);

        figure(42)
        clf
        hold on
        r.data.N1_time(r.data.time_key(counter)) = r.data.Rsum(i1,1);
        r.data.N2_time(r.data.time_key(counter)) = r.data.Rsum(i1,2);
        C2_mean = mean(r.data.N2_time,2,"omitnan");
        C2_std = std(r.data.N2_time,0,2,"omitnan");
        errorbar(x,C2_mean,C2_std,"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
        xlim([min(x) max(x)])
        ylim([0 1])
        xlabel('Pulse Time (microsec)')
        ylabel('N_2')

    elseif  r.c(1) > r.data.no_averages*(r.data.delta_param + r.data.time_param)

        counter = i1 - r.data.no_averages*(r.data.delta_param + r.data.time_param);

        % plotting as you go
        A1=r.data.phase_key';
        B1=r.data.Rsum([r.data.no_averages*(r.data.delta_param + r.data.time_param):i1],2);
        [A1,index(1:counter)] = sortrows(A1(1:counter));
        B1 = B1(index(counter),:);
        A2=r.data.phase_key';
        B2=r.data.Rsum([r.data.no_averages*(r.data.delta_param + r.data.time_param):i1],1);
        [A2,index(1:counter)] = sortrows(A2(1:counter));
        B2 = B2(index(counter),:);
        x = linspace(min(r.data.phase),max(r.data.phase),r.data.phase_param);

        figure(43)
        clf
        hold on
        r.data.N1_phase(r.data.phase_key(counter)) = r.data.Rsum(i1,1);
        r.data.N2_phase(r.data.phase_key(counter)) = r.data.Rsum(i1,2);
        C1_mean = mean(r.data.N1_phase,2,"omitnan");
        C1_std = std(r.data.N1_phase,0,2,"omitnan");
        C2_mean = mean(r.data.N2_phase,2,"omitnan");
        C2_std = std(r.data.N2_phase,0,2,"omitnan");
        errorbar(x,C1_mean - C2_mean,sqrt(C2_std.^2+C1_std.^2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
        xlim([min(x) max(x)])
        ylim([0 1])
        xlabel('Pulse Time (microsec)')
        ylabel('N_2')
    end

    if r.c(1) == r.data.no_averages*r.data.delta_param
        data = r.data;
        saveas(gcf,sprintf('D:\\data\\VMG_autosave\\Two Photon.fig'));
        save('D:\data\VMG_autosave\two_photon_data.mat','data');
        r.data.two_photon_opt = x(find(C2_mean == max(C2_mean)));
% r.data.two_photon_opt = 20.025;
        r.data.counter = i1;
    elseif r.c(1) == r.data.no_averages*(r.data.delta_param + r.data.time_param)
        data = r.data;
        saveas(gcf,sprintf('D:\\data\\VMG_autosave\\Time.fig'));
        save('D:\data\VMG_autosave\time_data.mat','data');
        r.data.time_opt = x(find(C2_mean == max(C2_mean)));
        r.data.counter = i1;
    elseif r.c(1) == r.data.no_averages*(r.data.delta_param + r.data.time_param + r.data.phase_param)
%         [f,g] = createFit(x',(0.5*(1 + C1_mean - C2_mean))');
%         x=linspace(0,360);
%         y = f.a*sin(pi*x/180+f.b)+f.c;
%         plot(x,y,'Linewidth',2,'Color',"blue");
%         r.data.vis = f.a;
        data = r.data;
        saveas(gcf,sprintf('D:\\data\\VMG_autosave\\Phase.fig'));
        save('D:\data\VMG_autosave\phase_data.mat','data');
    end

end