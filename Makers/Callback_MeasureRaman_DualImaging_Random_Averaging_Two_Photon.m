function Callback_MeasureRaman_DualImaging_Random_Averaging_Two_Photon(r)

if r.isInit()

    r.data.no_averages=10;
    r.data.number_of_param=26; %101
    r.data.abs_analysis = 1; %switches Abs_Analysis off (0) or on (1)

    r.data.param = repmat(linspace(19.9,20.15,r.data.number_of_param),1,r.data.no_averages); %freq
    %         r.data.param = repmat(linspace(0,pi,r.data.number_of_param),1,r.data.no_averages); %SLM phase
    %                 r.data.param = repmat(linspace(0,180,r.data.number_of_param),1,r.data.no_averages); %AOM phase
    %     r.data.param = repmat(linspace(2,100,r.data.number_of_param)*1e-6,1,r.data.no_averages); %pulse time
    %  r.data.param = repmat(linspace(0,1,r.data.number_of_param),1,r.data.no_averages); %mag bias
    %      r.data.param = repmat(linspace(0,1,r.data.number_of_param),1,r.data.no_averages); %power

    r.data.key=randperm(r.data.number_of_param*r.data.no_averages); %getting random key
    r.data.param = r.data.param(r.data.key); %randomising using key
    r.data.limit = 10; %error counter

    r.c.setup('var',r.data.param);
elseif r.isSet()
    r.make(r.devices.opt,'params','params',r.data.param(r.c(1))).upload;
    %     fprintf(1,'Run %d/%d, Time = %.3f µs\n',r.c.now,r.c.total,1e6*r.data.param(r.c(1)));
    %      fprintf(1,'Run %d/%d, Detuning = %.3f MHz\n',r.c.now,r.c.total,r.data.param(r.c(1)));
    fprintf(1,'Run %d/%d, delta = %.3f MHz \n',r.c.now,r.c.total,r.data.param(r.c(1)));

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

        %                 if (~img(1).raw.status.ok() || r.data.N(i1,2) == 0)
        %                             if (~img(1).raw.status.ok() || r.data.N(i1,1) == 0)
        %
        %                     Checks for an error in loading the files (caused by a missed
        %                     image) and reruns the last sequence
        %
        %                     if (~img(1).raw.status.ok())
        %                         warning('Imaging failed!')
        %                         r.c.decrement;
        %
        %                     elseif r.data.N(i1,2) == 0
        %                                     elseif (r.data.N(i1,1) == 0)
        %                         r.data.limit=1+r.data.limit
        %                         warning('Number of atoms too low!')
        %                         r.c.decrement;
        %                         if r.data.limit == 20
        %                             r.data.limit = 0;
        %                             r.stop
        %                             error('Run failed! lock your lasers.')
        %                         end
        %                     end
        %                     return;
        %                 elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        %                     r.c.decrement;
        %                     return;
        %                 end

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

    % plotting fringes as you go
    A1=r.data.key';
    B1=r.data.Rsum(:,2);
    [A1,index(1:r.c.i)] = sortrows(A1(r.c.i));
    B1 = B1(index(r.c.i),:);
    A2=r.data.key';
    B2=r.data.Rsum(:,1);
    [A2,index(1:r.c.i)] = sortrows(A2(1:r.c.i));
    B2 = B2(index(r.c.i),:);

    if r.c.i == 1
        r.data.D1 = NaN*ones(r.data.number_of_param,r.data.no_averages);
        r.data.D2 = r.data.D1;
        r.data.x = linspace(min(r.data.param),max(r.data.param),r.data.number_of_param);
    end

    %for phase
    figure(42)
    clf
    hold on
    r.data.D1(r.data.key(r.c.i)) = r.data.Rsum(i1,1);
    r.data.D2(r.data.key(r.c.i)) = r.data.Rsum(i1,2);
    C1_mean = mean(r.data.D1,2,"omitnan");
    C1_std = std(r.data.D1,0,2,"omitnan");
    C2_mean = mean(r.data.D2,2,"omitnan");
    C2_std = std(r.data.D2,0,2,"omitnan");
    y1=[(C1_mean-C2_mean)];
    %     y1=[(C2_mean)];
    %         errorbar(r.data.x,y1,sqrt(C1_std.^2+C2_std.^2)./sum(1-isnan(r.data.D1),2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    errorbar(r.data.x,0.5*(1+y1),sqrt(C1_std.^2+C2_std.^2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    %     errorbar(r.data.x,y1,C2_std,"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    xlim([min(r.data.x) max(r.data.x)])
    %         xticks([0 45 90 135 180])
    ylim([0 1])
    xlabel('Two Photon Detuning (MHz)')
    ylabel('N_2')
    %             if r.c(1) > 20
    %                 [f,g] = createFit(r.data.x', y1');
    %                 x=linspace(0,180);
    %                 y = f.a*sin(2*pi*x/180+f.b)+f.c;
    %                 plot(x,y,'Linewidth',2,'Color',"blue");
    %             end

    % %for time
    %     figure(42)
    %     clf
    %     hold on
    %     r.data.D1(r.data.key(r.c.i)) = r.data.Rsum(i1,1);
    %     r.data.D2(r.data.key(r.c.i)) = r.data.Rsum(i1,2);
    %     C1_mean = mean(r.data.D1,2,"omitnan");
    %     C1_std = std(r.data.D1,0,2,"omitnan");
    %     C2_mean = mean(r.data.D2,2,"omitnan");
    %     C2_std = std(r.data.D2,0,2,"omitnan");
    %     y1=C2_mean*100;
    %     %         errorbar(r.data.x*1e6,y1,C2_std./sum(1-isnan(r.data.D2),2),"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    %     errorbar(r.data.x*1e6,y1,C2_std*100,"o","MarkerSize",5,"MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90],'linewidth',2,'Color','blue')
    %     xlim([min(r.data.x)*1e6 max(r.data.x)*1e6])
    %     ylim([0 50])
    %     xlabel('Pulse Time (\mus)')
    %     ylabel('N_2 Transfer')
    % %     if r.c(1) > 10
    % %         [f,g] = createFit_time(r.data.x', y1');
    % %         x=linspace(0,50,51)*1e-6;
    % %         y = f.a*exp(f.b*x).*sin(f.c*x+f.d)+f.g;
    % %         plot(x*1e6,y,'Linewidth',2,'Color',"blue");
    % %     end

    if r.c.done(1)

        %saving to VMG_autosave folder in D drive
        data = r.data;
        save('D:\data\VMG_autosave\two_photon_data.mat','data');
        saveas(gcf,'D:\data\VMG_autosave\two_photon_data.fig')

    end
end


end