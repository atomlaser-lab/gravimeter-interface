function Callback_RH(r)
%% Inuts
Param = 8:2:10;
NumAverage = 2;

if r.isInit()
    r.data.no_averages = NumAverage;
    r.data.number_of_param = numel(Param);

    r.data.param = repmat(Param,1,r.data.no_averages); %freq
    r.data.key = randperm(numel(Param)*r.data.no_averages); %getting random key
    r.data.param = r.data.param(r.data.key); %randomising using key
    r.data.limit = 1; %error counter
    r.c.setup('var',r.data.param);



elseif r.isSet()
    r.make(r.devices.opt,'params','params',r.data.param(r.c(1))).upload;
    fprintf(1,'Run %d/%d, Rotation = %.3f degres \n',r.c.now,r.c.total,r.data.param(r.c(1)));
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1 + 0.5*rand);
    if r.data.abs_analysis==1
        img = Abs_Analysis_DualState('last',1);
        %extract data
        r.data.files{i1,1} = img(1).raw.files;
        r.data.N(i1,:) = img.get('N');
        r.data.Nsum(i1,:) = img.get('Nsum');
        r.data.peakOD(i1,:) = img.get('peakOD');
        r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:),2);
        r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:),2);

        if (~img(1).raw.status.ok())
            %
            % Checks for an error in loading the files (caused by a missed
            % image) and reruns the last sequence
            %
            if (~img(1).raw.status.ok())
                warning('Imaging failed!')
                r.c.decrement;
            end
            return;
        elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
            r.c.decrement;
            return;
        end
    end
    % sort param
    A1=r.data.key';
    B1=r.data.Rsum(:,2);
    [A1,index(1:r.c.i)] = sortrows(A1(r.c.i));
%     sort data to param
    B1 = B1(index(r.c.i),:);
    A2=r.data.key';
    B2=r.data.Rsum(:,1);
    [A2,index(1:r.c.i)] = sortrows(A2(1:r.c.i));
    B2 = B2(index(r.c.i),:);
    if r.c.i == 1
        % create nan vector for all possible runs
        r.data.D1 = NaN*ones(r.data.number_of_param,r.data.no_averages);
        r.data.D2 = r.data.D1;
        r.data.x = linspace(min(r.data.param),max(r.data.param),r.data.number_of_param);
    end
    %for phase
    figure(42)
    clf
    hold on
%     populate nan vector with sorted vector
    r.data.D1(r.data.key(r.c.i)) = r.data.Rsum(i1,1);
    r.data.D2(r.data.key(r.c.i)) = r.data.Rsum(i1,2);
    % calculate the mean for each variable
    C1_mean = mean(r.data.D1,2,"omitnan");
    C1_std = std(r.data.D1,0,2,"omitnan");
    C2_mean = mean(r.data.D2,2,"omitnan");
    C2_std = std(r.data.D2,0,2,"omitnan");
    %         y1=[(C1_mean-C2_mean)];
    y1=[(C2_mean)];
    %         errorbar(r.data.x,y1,sqrt(C1_std.^2+C2_std.^2)./sum(1-isnan(r.data.D1),2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    errorbar(r.data.x,y1,sqrt(C1_std.^2+C2_std.^2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    xlim([min(r.data.x) max(r.data.x)])
    %         xticks([0 45 90 135 180])
    ylim([0 1])
    xlabel('LG Phase Rotation')
    ylabel('N_2-N_1')

    if r.c.done(1)
        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\data.mat','data');
        if r.data.abs_analysis==1
            %
            %sorting data
            A1=data.key';
            B1=data.Rsum(:,1);
            [A1,index] = sortrows(A1);
            B1 = B1(index,:);
            B1 = reshape(B1,data.number_of_param,data.no_averages);
            %             B1(16,1) = 0.35;
            C1_mean = mean(B1,2);
            C1_std = std(B1,0,2);
            A2=data.key';
            B2=data.Rsum(:,2);
            [A2,index] = sortrows(A2);
            B2 = B2(index,:);
            B2 = reshape(B2,data.number_of_param,data.no_averages);
            %             B1(16,1) = 0.35;
            C2_mean = mean(B2,2);
            C2_std = std(B2,0,2);

            clf
            hold on
            set(gcf,'color','w');
            x=linspace(min(data.param),max(data.param),data.number_of_param);
            y1=[100-C1_mean*100];
            errorbar(x*1.6*23.5,y1,[C1_std]*100,"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')

        end
    end
end