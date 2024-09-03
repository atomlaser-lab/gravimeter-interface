function Callback_Rhys_2State(r)

% % Inputs
FigNum = 5;

TitleStuff.TotalPower = '20';
TitleStuff.P_rat = '7';
TitleStuff.Mag = '1.5';
TitleStuff.t_0 = '10';
TitleStuff.T = '0.5';
TitleStuff.Tau  = '5';
TitleStuff.SPD = '4.95';
TitleStuff.TPD = '-20 + 33*1e-3';

TitleStuff.SubTitle = 'UD = 10 A, EW = 20 V, 3-Pulse';
Title = append('P_{total} = ',TitleStuff.TotalPower,' mW, ','P_S/P_C = ',TitleStuff.P_rat,', ', TitleStuff.Mag,'x Mag, ', 't_0 = ',TitleStuff.t_0,' us, ', '\Delta = ',TitleStuff.SPD,' GHz,','\delta = ',TitleStuff.TPD,' MHz, ','\tau = ',TitleStuff.Tau,' us', ', T = ',TitleStuff.T,' ms');
TitleStuff.Title = Title;

% Title = 'In-Trap MW';
% TitleStuff.SubTitle = 'pi stability of pulse 1';

% Param = const.randomize(unique([0:10:360]));
Param = [10:10:200];
PlotParam = 1*Param;
ParamName = ScanableParameters.TwoPhoton;
% ParamName = 'I_{sideband}/I_{carrier}';

% % % If there are multiple ROIs, what do you want to count?
F2_ROI = 3;
F1_ROI = 1;
allROI = 0;
MovAv = 0;


if r.isInit()
    r.data.Param = Param;
    r.data.PlotParam = PlotParam;
    r.c.setup('var',r.data.Param);
elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.Param(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, T = %.0f us\n',r.c.now,r.c.total,...
        r.data.Param(r.c(1)));

elseif r.isAnalyze()
    i1 = r.c(1);
    %     pause(0.5 + 0.5*rand);
    img = Abs_Analysis_DualState_RT('last');

    if r.c(1) == 1
        % Create a structure for each ROI for each image
        F2field_names = arrayfun(@(x) sprintf('ROI%d', x), 1:numel(img(1).clouds), 'UniformOutput', false);
        F1field_names = arrayfun(@(x) sprintf('ROI%d', x), 1:numel(img(2).clouds), 'UniformOutput', false);
        r.data.F2 = cell2struct(cell(size(F2field_names)), F2field_names, 2);
        r.data.F1 = cell2struct(cell(size(F1field_names)), F1field_names, 2);
    end

    if ~img(1).raw.status.ok()

        %         Checks for an error in loading the files (caused by a missed
        %         image) and reruns the last sequence

        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end

    %     Store raw data

    r.data.files{i1,1} = img(1).raw.files;

    %
    % % % ROI error Check
    %

    if F2_ROI > size(img(1).clouds,1)
        F2_ROI = size(img(1).clouds,1);
    end
    if F1_ROI > size(img(2).clouds,1)
        F1_ROI = size(img(2).clouds,1);
    end

    %
    % Get processed data for input region of interest
    %
    r.data.N(i1,:) = [img(1).clouds(F2_ROI).N,img(2).clouds(F1_ROI).N];
    r.data.Nsum(i1,:) = [img(1).clouds(F2_ROI).Nsum,img(2).clouds(F1_ROI).Nsum];

    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));

        r.data.T_F1(i1,:) = [img(2).clouds(F1_ROI).T];
        r.data.T_F2(i1,:) = [img(1).clouds(F2_ROI).T];

        %
        % Get processed data for all regions of interest
        %
        if allROI == 1
            F2Name = fieldnames(r.data.F2);
            for ii = 1:size(img(1).clouds,1)
                r.data.F2.(F2Name{ii}).N(i1,:) = img(1).clouds(ii).N;
                r.data.F2.(F2Name{ii}).Nsum(i1,:) = img(1).clouds(ii).Nsum;
                r.data.F2.(F2Name{ii}).T(i1,:) = img(1).clouds(ii).T;
                r.data.F2.(F2Name{ii}).OD(i1) = img(1).clouds(ii).peakOD;
    
                r.data.F2.(F2Name{ii}).R(i1,:) = img(1).clouds(ii).N/(sum(vertcat(img(1).clouds(:).N)) + sum(vertcat(img(2).clouds(:).N)));
                r.data.F2.(F2Name{ii}).Rsum(i1,:) = img(1).clouds(ii).Nsum/(sum(vertcat(img(1).clouds(:).Nsum)) + sum(vertcat(img(2).clouds(:).Nsum)));
            end
    
            F1Name = fieldnames(r.data.F1);
            for ii = 1:size(img(2).clouds,1)
                r.data.F1.(F1Name{ii}).N(i1,:) = img(2).clouds(ii).N;
                r.data.F1.(F1Name{ii}).Nsum(i1,:) = img(2).clouds(ii).Nsum;
                r.data.F1.(F1Name{ii}).T(i1,:) = img(2).clouds(ii).T;
                r.data.F1.(F2Name{ii}).OD(i1) = img(2).clouds(ii).peakOD;
    
    
                r.data.F1.(F2Name{ii}).R(i1,:) = img(2).clouds(ii).N/(sum(vertcat(img(1).clouds(:).N)) + sum(vertcat(img(2).clouds(:).N)));
                r.data.F1.(F2Name{ii}).Rsum(i1,:) = img(2).clouds(ii).Nsum/(sum(vertcat(img(1).clouds(:).Nsum)) + sum(vertcat(img(2).clouds(:).Nsum)));
            end
        end
        
        % % % Grab Beam position IF image exists
        %
        if size(img(1).raw.images,3) == 6
            [Beam1_cy,Beam1_cx] = find_beam_position(img(1).raw.images(:,:,5));
            [Beam2_cy,Beam2_cx] = find_beam_position(img(1).raw.images(:,:,6));
    
            r.data.Beam1.y.amp(i1) = Beam1_cy(1);
            r.data.Beam1.y.pos(i1) = Beam1_cy(2);
            r.data.Beam1.y.width(i1) = Beam1_cy(3);
            r.data.Beam1.x.amp(i1) = Beam1_cx(1);
            r.data.Beam1.x.pos(i1) = Beam1_cx(2);
            r.data.Beam1.x.width(i1) = Beam1_cx(3);
    
            r.data.Beam2.y.amp(i1) = Beam2_cy(1);
            r.data.Beam2.y.pos(i1) = Beam2_cy(2);
            r.data.Beam2.y.width(i1) = Beam2_cy(3);
            r.data.Beam2.x.amp(i1) = Beam2_cx(1);
            r.data.Beam2.x.pos(i1) = Beam2_cx(2);
            r.data.Beam2.x.width(i1) = Beam2_cx(3);
        end


    %% Plots
    figure(FigNum);clf
    subplot(4,2,[1 3 5])
    scatter(r.data.PlotParam(1:i1),r.data.N(1:i1,:),'filled');
    plot_format(ParamName,'N_{F = i}','',12);
    grid on
    ylim([0,Inf]);
    legend('F2','F1')

    subplot(4,2,7)
    scatter(r.data.PlotParam(1:i1),r.data.N(1:i1,1) + r.data.N(1:i1,2),'filled');
    plot_format(ParamName,'N_{Total}','',12);

    subplot(4,2,[2 4])
    scatter(r.data.PlotParam(1:i1),r.data.Rsum(1:i1,1),100,'b');
    hold on
    scatter(r.data.PlotParam(1:i1),r.data.R(1:i1,1),40,'b','filled');
    plot_format(ParamName,'F = 2 Population','',12);
    grid on;

    subplot(4,2,[6 8])
    scatter(r.data.PlotParam(1:i1),r.data.Rsum(1:i1,2),100,'r');
    hold on
    scatter(r.data.PlotParam(1:i1),r.data.R(1:i1,2),40,'r','filled');
    plot_format(ParamName,'F = 1 Population','',12);
    grid on;

    if MovAv == 1
        subplot(4,2,[2 4])
        scatter(r.data.PlotParam(1:i1),movmean(r.data.Rsum(1:i1,1),5),20,'k','filled');
        legend('Sum','Fit','Moving Average')
        subplot(4,2,[6 8])
        scatter(r.data.PlotParam(1:i1),movmean(r.data.Rsum(1:i1,2),5),20,'k','filled');
    end



    sgtitle({['{\bf\fontsize{14}' Title '}'],TitleStuff.SubTitle});

    % % % % % % % % % % % %
    if allROI == 1
        maxlength = max(numel(img(1).clouds),numel(img(2).clouds));
        figure(FigNum+1);clf
        sgtitle({['{\bf\fontsize{14}' Title '}'],TitleStuff.SubTitle});
        for ii = 1:numel(img(1).clouds)
            subplot(maxlength,2,ii*2-1);
            if i1 == 1
                hold on;
            end
            scatter(r.data.PlotParam(1:i1), r.data.F2.(F2Name{ii}).Rsum,'filled')
            hold on
            scatter(r.data.PlotParam(1:i1),r.data.F2.(F2Name{ii}).R,'filled')

            ylabel('Pop')
            if ii == 1
                title(sprintf('F = 2 atoms \n ROI %g',ii))
            else
                title(sprintf('ROI %g',ii))
            end
            if ii == numel(img(1).clouds)
                xlabel(ParamName)
            end
        end
        for ii = 1:numel(img(2).clouds)
            subplot(maxlength,2,ii*2);
            if i1 == 1
                hold on;
            end
            scatter(r.data.PlotParam(1:i1),r.data.F1.(F1Name{ii}).Rsum,'filled')
            hold on
            scatter(r.data.PlotParam(1:i1),r.data.F1.(F1Name{ii}).R,'filled')

            if ii == 1
                title(sprintf('F = 1 atoms \n ROI %g',ii))
            else
                title(sprintf('ROI %g',ii))
            end
            ylabel('Pop')
            if ii == numel(img(1).clouds)
                xlabel(ParamName)
            end
        end
    end

    % % % % % % % % % % % % % %
    if size(img(1).raw.images,3) == 6

        figure(FigNum + 2);clf
        % Beam 1
        subplot(3,2,1)
        scatter(r.data.PlotParam(1:i1),r.data.Beam1.y.amp,'r','filled')
        hold on
        scatter(r.data.PlotParam(1:i1),r.data.Beam2.y.amp,'b','filled')
        ylabel('amp')
        title('Vertical')
        legend('Beam 1','Beam 2')

        subplot(3,2,3)
        scatter(r.data.PlotParam(1:i1),r.data.Beam1.y.pos,'r','filled')
        hold on
        scatter(r.data.PlotParam(1:i1),r.data.Beam2.y.pos,'b','filled')
        ylabel('pos')

        subplot(3,2,5)
        scatter(r.data.PlotParam(1:i1),r.data.Beam1.y.width,'r','filled')
        hold on
        scatter(r.data.PlotParam(1:i1),r.data.Beam2.y.width,'b','filled')
        ylabel('width')
        xlabel(ParamName)


        % Beam 2
        subplot(3,2,2)
        scatter(r.data.PlotParam(1:i1),r.data.Beam1.x.amp,'r','filled')
        hold on
        scatter(r.data.PlotParam(1:i1),r.data.Beam2.x.amp,'b','filled')
        ylabel('amp')
        title('Horizontal')
        legend('Beam 1','Beam 2')

        subplot(3,2,4)
        scatter(r.data.PlotParam(1:i1),r.data.Beam1.x.pos,'r','filled')
        hold on
        scatter(r.data.PlotParam(1:i1),r.data.Beam2.x.pos,'b','filled')
        ylabel('pos')

        subplot(3,2,6)
        scatter(r.data.PlotParam(1:i1),r.data.Beam1.x.width,'r','filled')
        hold on
        scatter(r.data.PlotParam(1:i1),r.data.Beam2.x.width,'b','filled')
        ylabel('width')
        xlabel(ParamName)
    end
end
