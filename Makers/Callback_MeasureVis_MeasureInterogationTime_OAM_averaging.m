function Callback_MeasureVis_MeasureInterogationTime_OAM_averaging(r)

if r.isInit()
    r.data.T = [1e-6 10e-6 100e-6 1e-3]; %for dipoles: [1.5:0.5:8]
    r.data.OAM = [1 2 5 10];
    r.data.no_averages = 5;
    r.data.number_of_param = 5;
    r.data.phase = repmat(linspace(0,180,r.data.number_of_param),1,r.data.no_averages); % phase AOM
    r.data.limit = 1;
    r.data.vis = NaN*ones(length(r.data.OAM),length(r.data.T));
    r.data.key = randperm(r.data.number_of_param*r.data.no_averages); %getting random key
    r.data.phase = r.data.phase(r.data.key);
    r.c.setup('var',r.data.phase,r.data.T,r.data.OAM);
elseif r.isSet()
    r.make(r.devices.opt,'param1',r.data.phase(r.c(1)),'param2',r.data.T(r.c(2)),'param3',r.data.OAM(r.c(3))).upload;
    fprintf(1,'Run %d/%d, Rotation = %.3f degres, T = %.3f μs, OAM = %d\n',r.c.now,r.c.total,2*r.data.phase(r.c(1)),r.data.T(r.c(2))*1e6,r.data.OAM(r.c(3)));
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    i3 = r.c(3);
    pause(0.1 + 0.5*rand);
    img = Abs_Analysis_DualState('last',1);

    r.data.files{i1,i2,i3,1} = img(1).raw.files;
    r.data.N(i1,i2,i3,:) = img.get('N');
    r.data.Nsum(i1,i2,i3,:) = img.get('Nsum');
    r.data.peakOD(i1,i2,i3,:) = img.get('peakOD');
    r.data.R(i1,i2,i3,:) = r.data.N(i1,i2,i3,:)./sum(r.data.N(i1,i2,i3,:));
    r.data.Rsum(i1,i2,i3,:) = r.data.Nsum(i1,i2,i3,:)./sum(r.data.Nsum(i1,i2,i3,:));

    if (~img(1).raw.status.ok() || r.data.N(i1,i2,i3,1) == 0)
        %         if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0)
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        if (~img(1).raw.status.ok())
            warning('Imaging failed!')
            r.c.decrement;

        elseif r.data.N(i1,i2,i3,1) == 0
            %             elseif (r.data.N(i1,1) == 0)
            r.data.limit=1+r.data.limit;
            warning('Number of atoms too low!')
            r.c.decrement;
            if r.data.limit == 20
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

    % plotting fringes as you go
    A1=r.data.key';
    B1=r.data.Rsum(:,i2,i3,2);
    [A1,index(1:i1)] = sortrows(A1(1:i1));
    B1 = B1(index(i1),:);
    A2=r.data.key';
    B2=r.data.Rsum(:,i2,i3,1);
    [A2,index(1:i1)] = sortrows(A2(1:i1));
    B2 = B2(index(i1),:);

    if i1 == 1
        r.data.D1 = NaN*ones(r.data.number_of_param,r.data.no_averages);
        r.data.D2 = r.data.D1;
        r.data.x = 2*linspace(min(r.data.phase),max(r.data.phase),r.data.number_of_param);
    end

    figure(42)
    clf
    hold on
    r.data.D1(r.data.key(i1)) = r.data.Rsum(i1,i2,i3,1);
    r.data.D2(r.data.key(i1)) = r.data.Rsum(i1,i2,i3,2);
    C1_mean = mean(r.data.D1,2,"omitnan");
    C1_std = std(r.data.D1,0,2,"omitnan");
    C2_mean = mean(r.data.D2,2,"omitnan");
    C2_std = std(r.data.D2,0,2,"omitnan");
    y1=[(C1_mean-C2_mean)];
    errorbar(r.data.x,0.5*(1+y1),sqrt(C1_std.^2+C2_std.^2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    xlim([min(r.data.x) max(r.data.x)])
    ylim([0 1])
    xlabel('LG Rotation Phase (Deg)')
    ylabel('N_2-N1')

%     figure(133);clf;
%     hold on
%     x = 2*r.data.phase;
%     if i1 == 1
%         r.data.y2 = NaN*ones(1,length(x));
%     end
%     r.data.y2(i1) = 0.5*(1+r.data.Rsum(i1,i2,1)-r.data.Rsum(i1,i2,2));
%     plot(x,r.data.y2,'o');
%     enhformat('Rotated Phase','Population Difference')
%     grid on;
%     ylim([0 1])
%     xlim([0 360])

    if r.c(1) == length(r.data.phase)
        [f,g] = createFit(r.data.x',(0.5*(1+y1))');
        x=linspace(0,360);
        y = f.a*sin(pi*x/180+f.b)+f.c;
        plot(x,y,'Linewidth',2,'Color',"blue");
        saveas(gcf,sprintf('D:\\data\\VMG_autosave\\T = %.1f μs, OAM = %d.png',r.data.T(i2)*1e6,r.data.OAM(i3)));
        r.data.vis(i3,i2) = 2*f.a;

        figure(43)
        clf
        hold on
        plot(r.data.T*1e6,r.data.vis(1,:),'o')
        plot(r.data.T*1e6,r.data.vis(2,:),'o')
        plot(r.data.T*1e6,r.data.vis(3,:),'o')
        plot(r.data.T*1e6,r.data.vis(4,:),'o')
        ylim([0 max(max(r.data.vis))])
        grid on;
        enhformat('Interrogation Time (μs)','Visibility')
    end

    if r.c.done(1)

        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\data.mat','data');
        saveas(gcf,'D:\data\VMG_autosave\data.fig');

    end

end