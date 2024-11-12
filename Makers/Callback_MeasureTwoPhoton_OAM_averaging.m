function Callback_MeasureTwoPhoton_OAM_averaging(r)

if r.isInit()
    r.data.OAM = [1e-6 10e-6 100e-6 200e-6 500e-6]; %[0 1 2 5 10] %for dipoles: [1.5:0.5:8] %for RF: [0.5:0.25:2]
    r.data.no_averages = 1;
    r.data.number_of_param = 20;
    r.data.delta = repmat(linspace(0,180,r.data.number_of_param),1,r.data.no_averages); % phase AOM
    r.data.limit = 100;
    r.data.key = randperm(r.data.number_of_param*r.data.no_averages); %getting random key
    r.data.delta = r.data.delta(r.data.key);
    r.c.setup('var',r.data.delta,r.data.OAM);
elseif r.isSet()
    r.make(r.devices.opt,'param1',r.data.delta(r.c(1)),'param2',r.data.OAM(r.c(2))).upload;
    fprintf(1,'Run %d/%d, LG Rotation = %.3f degrees, OAM = %d\n',r.c.now,r.c.total,r.data.delta(r.c(1)),r.data.OAM(r.c(2)));
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

     if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0 || r.data.N(i1,1) > 1e7  || (r.data.N(i1,1) == 0 && r.data.N(i1,2) == 0))
      %         if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0)
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        if (~img(1).raw.status.ok())
            warning('Imaging failed!')
            r.c.decrement;

        elseif (r.data.N(i1,1) == 0 || r.data.N(i1,1) > 1e7 || (r.data.N(i1,1) == 0 && r.data.N(i1,2) == 0))
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
    B1=r.data.Rsum(:,i2,2);
    [A1,index(1:i1)] = sortrows(A1(1:i1));
    B1 = B1(index(i1),:);
    A2=r.data.key';
    B2=r.data.Rsum(:,i2,1);
    [A2,index(1:i1)] = sortrows(A2(1:i1));
    B2 = B2(index(i1),:);

    if i1 == 1
        r.data.D1 = NaN*ones(r.data.number_of_param,r.data.no_averages);
        r.data.D2 = r.data.D1;
        r.data.x = linspace(min(r.data.delta),max(r.data.delta),r.data.number_of_param);
    end

    figure(42)
    clf
    hold on
    r.data.D1(r.data.key(i1)) = r.data.Rsum(i1,i2,1);
    r.data.D2(r.data.key(i1)) = r.data.Rsum(i1,i2,2);
    C1_mean = mean(r.data.D1,2,"omitnan");
    C1_std = std(r.data.D1,0,2,"omitnan");
    C2_mean = mean(r.data.D2,2,"omitnan");
    C2_std = std(r.data.D2,0,2,"omitnan");
    y1=[(C2_mean - C1_mean)];
%  y1=[C2_mean];
    scatter(r.data.x,y1,'filled')
%     errorbar(r.data.x,y1,sqrt(C2_std.^2+C1_std.^2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    xlim([min(r.data.x) max(r.data.x)])
    ylim([-1 1])
    xlabel('LG Rotation (deg)')
    ylabel('Population Difference')

    if r.c(1) == length(r.data.delta)
          saveas(gcf,sprintf('D:\\data\\VMG_autosave\\VMG for OAM =  %d.fig',r.data.OAM(i2)));
    end

    if r.c.done(1)

        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\data.mat','data');
        saveas(gcf,'D:\data\VMG_autosave\data.fig');

    end

end