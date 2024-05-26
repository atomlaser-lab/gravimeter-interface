function Callback_Rhys_2State_Average(r)

% % Inputs
ClearImage = 1;
FigNum = 5;

Title = 'Pumping: P_{total} = 1 mW, P_S/P_C = 7/1, 3x Mag, t_0 = 0 us, \Delta = 4.95 GHz, \delta = -20 + param, \tau = 240 us';
SubTitle = 'EW = scan, NS = 0, UD = 0';
NumAverages = 7;
Param = -.75:0.25:.75;


PlotFactor = 1;
ParamName = ScanableParameters.TwoPhoton;
% ParamName = 'EW Bias (V)';
% % % If there are multiple ROIs, what do you want to count?
F2_ROI = 3;
F1_ROI = 2;

% % Sequence
if r.isInit()
    % % % % Randomise data
    r.data.RandomOrder = randperm(numel(Param)*NumAverages);
    RepeatedParam = repmat(Param,1,NumAverages);
    RandomData = RepeatedParam(r.data.RandomOrder);

    % % % Key to Sort Data
    [r.data.ParamIndex, r.data.AvIndex] = ind2sub([numel(Param),NumAverages],r.data.RandomOrder);

    % % % % Stored data
    r.data.Param = RandomData;
    r.data.PlotParam = RandomData*PlotFactor;
    r.c.setup('var',r.data.Param);
    if ClearImage == 1
        figure(FigNum);clf
        figure(FigNum+1);clf
    end

elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.Param(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, T = %.0f us\n',r.c.now,r.c.total,...
        r.data.Param(r.c(1)));

elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis_DualState_RT('last');

    if r.c(1) == 1
        % Create a structure for each ROI for each image
        F2field_names = arrayfun(@(x) sprintf('ROI%d', x), 1:numel(img(1).clouds), 'UniformOutput', false);
        F1field_names = arrayfun(@(x) sprintf('ROI%d', x), 1:numel(img(2).clouds), 'UniformOutput', false);
        r.data.F2 = cell2struct(cell(size(F2field_names)), F2field_names, 2);
        r.data.F1 = cell2struct(cell(size(F1field_names)), F1field_names, 2);
    end

    %     if ~img(1).raw.status.ok()
    %
    % %         Checks for an error in loading the files (caused by a missed
    % %         image) and reruns the last sequence
    %
    %         r.c.decrement;
    %         return;
    %     elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
    %         r.c.decrement;
    %         return;
    %     end

    % Store raw data
    r.data.Unsorted.files{i1,1} = img(1).raw.files;

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
    AvIndex = r.data.AvIndex(r.c.i);
    ParIndex = r.data.ParamIndex(r.c.i);


    if r.c.i == 1
        % Pre-allocate data
        r.data.All.N = nan(numel(Param),2,NumAverages);
        r.data.All.R = nan(numel(Param),2,NumAverages);
        r.data.All.NSum = nan(numel(Param),2,NumAverages);
        r.data.All.RSum = nan(numel(Param),2,NumAverages);
    end

    % All data
    r.data.All.N(ParIndex,:,AvIndex) = [img(1).clouds(F2_ROI).N,img(2).clouds(F1_ROI).N];
    r.data.All.NSum(ParIndex,:,AvIndex) = [img(1).clouds(F2_ROI).Nsum,img(2).clouds(F1_ROI).Nsum];
    r.data.All.R(ParIndex,:,AvIndex) = r.data.All.N(ParIndex,:,AvIndex)./sum(r.data.All.N(ParIndex,:,AvIndex));
    r.data.All.RSum(ParIndex,:,AvIndex) = r.data.All.NSum(ParIndex,:,AvIndex)./sum(r.data.All.NSum(ParIndex,:,AvIndex));
    % Averaged data
    r.data.N = nanmean(r.data.All.N,3);
    r.data.NSum = nanmean(r.data.All.NSum,3);
    r.data.R = nanmean(r.data.All.R,3);
    r.data.RSum = nanmean(r.data.All.RSum,3);
    % Std
    r.data.N_std = nanstd(r.data.All.N,0,3);
    r.data.NSum_std = nanstd(r.data.All.NSum,0,3);
    r.data.R_std = nanstd(r.data.All.R,0,3);
    r.data.RSum_std = nanstd(r.data.All.RSum,0,3);

    %
    % Get processed data for all regions of interest
    %
    F2Name = fieldnames(r.data.F2);
    for ii = 1:size(img(1).clouds,1)
        if r.c.i == 1
            % Create a nan matrix that gets populated
            r.data.F2.(F2Name{ii}).All.N = nan(numel(Param),NumAverages);
            r.data.F2.(F2Name{ii}).All.Nsum = nan(numel(Param),NumAverages);
            r.data.F2.(F2Name{ii}).All.T = nan(numel(Param),NumAverages,2);
            r.data.F2.(F2Name{ii}).All.OD = nan(numel(Param),NumAverages);

            r.data.F2.(F2Name{ii}).All.R = nan(numel(Param),NumAverages);
            r.data.F2.(F2Name{ii}).All.Rsum = nan(numel(Param),NumAverages);
        end

        % All data
        r.data.F2.(F2Name{ii}).All.N(ParIndex,AvIndex) = img(1).clouds(ii).N;
        r.data.F2.(F2Name{ii}).All.Nsum(ParIndex,AvIndex) = img(1).clouds(ii).Nsum;
        r.data.F2.(F2Name{ii}).All.T(ParIndex,AvIndex,:) = img(1).clouds(ii).T;
        r.data.F2.(F2Name{ii}).All.OD(ParIndex,AvIndex) = img(1).clouds(ii).peakOD;

        r.data.F2.(F2Name{ii}).All.R(ParIndex,AvIndex) = img(1).clouds(ii).N/(sum(vertcat(img(1).clouds(:).N)) + sum(vertcat(img(2).clouds(:).N)));
        r.data.F2.(F2Name{ii}).All.Rsum(ParIndex,AvIndex) = img(1).clouds(ii).Nsum/(sum(vertcat(img(1).clouds(:).Nsum)) + sum(vertcat(img(2).clouds(:).Nsum)));

        % Average Data
        r.data.F2.(F2Name{ii}).N = nanmean(r.data.F2.(F2Name{ii}).All.N,2);
        r.data.F2.(F2Name{ii}).Nsum = nanmean(r.data.F2.(F2Name{ii}).All.Nsum,2);
        r.data.F2.(F2Name{ii}).T = nanmean(r.data.F2.(F2Name{ii}).All.T,2);
        r.data.F2.(F2Name{ii}).OD = nanmean(r.data.F2.(F2Name{ii}).All.OD,2);

        r.data.F2.(F2Name{ii}).R = nanmean(r.data.F2.(F2Name{ii}).All.R,2);
        r.data.F2.(F2Name{ii}).Rsum = nanmean(r.data.F2.(F2Name{ii}).All.Rsum,2);
        % Std
        r.data.F2.(F2Name{ii}).N_std = nanstd(r.data.F2.(F2Name{ii}).All.N,0,2);
        r.data.F2.(F2Name{ii}).Nsum_std = nanstd(r.data.F2.(F2Name{ii}).All.Nsum,0,2);
        r.data.F2.(F2Name{ii}).T_std = nanstd(r.data.F2.(F2Name{ii}).All.T,0,2);
        r.data.F2.(F2Name{ii}).OD_std = nanstd(r.data.F2.(F2Name{ii}).All.OD,0,2);

        r.data.F2.(F2Name{ii}).R_std = nanstd(r.data.F2.(F2Name{ii}).All.R,0,2);
        r.data.F2.(F2Name{ii}).Rsum_std = nanstd(r.data.F2.(F2Name{ii}).All.Rsum,0,2);
    end

    F1Name = fieldnames(r.data.F1);
    for ii = 1:size(img(2).clouds,1)
        if r.c.i == 1
            % Create a nan matrix that gets populated
            r.data.F1.(F1Name{ii}).All.N = nan(numel(Param),NumAverages);
            r.data.F1.(F1Name{ii}).All.Nsum = nan(numel(Param),NumAverages);
            r.data.F1.(F1Name{ii}).All.T = nan(numel(Param),NumAverages,2);
            r.data.F1.(F1Name{ii}).All.OD = nan(numel(Param),NumAverages);

            r.data.F1.(F1Name{ii}).All.R = nan(numel(Param),NumAverages);
            r.data.F1.(F1Name{ii}).All.Rsum = nan(numel(Param),NumAverages);
        end
        % All Data
        r.data.F1.(F1Name{ii}).All.N(ParIndex,AvIndex) = img(2).clouds(ii).N;
        r.data.F1.(F1Name{ii}).All.Nsum(ParIndex,AvIndex) = img(2).clouds(ii).Nsum;
        r.data.F1.(F1Name{ii}).All.T(ParIndex,AvIndex,:) = img(2).clouds(ii).T;
        r.data.F1.(F1Name{ii}).All.OD(ParIndex,AvIndex) = img(2).clouds(ii).peakOD;

        r.data.F1.(F1Name{ii}).All.R(ParIndex,AvIndex,:) = img(2).clouds(ii).N/(sum(vertcat(img(1).clouds(:).N)) + sum(vertcat(img(2).clouds(:).N)));
        r.data.F1.(F1Name{ii}).All.Rsum(ParIndex,AvIndex,:) = img(2).clouds(ii).Nsum/(sum(vertcat(img(1).clouds(:).Nsum)) + sum(vertcat(img(2).clouds(:).Nsum)));

        % Average Data
        r.data.F1.(F1Name{ii}).N = nanmean(r.data.F1.(F1Name{ii}).All.N,2);
        r.data.F1.(F1Name{ii}).Nsum = nanmean(r.data.F1.(F1Name{ii}).All.Nsum,2);
        r.data.F1.(F1Name{ii}).T = nanmean(r.data.F1.(F1Name{ii}).All.T,2);
        r.data.F1.(F1Name{ii}).OD = nanmean(r.data.F1.(F1Name{ii}).All.OD,2);

        r.data.F1.(F1Name{ii}).R = nanmean(r.data.F1.(F1Name{ii}).All.R,2);
        r.data.F1.(F1Name{ii}).Rsum = nanmean(r.data.F1.(F1Name{ii}).All.Rsum,2);
        % Std
        r.data.F1.(F1Name{ii}).N_std = nanstd(r.data.F1.(F1Name{ii}).All.N,0,2);
        r.data.F1.(F1Name{ii}).Nsum_std = nanstd(r.data.F1.(F1Name{ii}).All.Nsum,0,2);
        r.data.F1.(F1Name{ii}).T_std = nanstd(r.data.F1.(F1Name{ii}).All.T,0,2);
        r.data.F1.(F1Name{ii}).OD_std = nanstd(r.data.F1.(F1Name{ii}).All.OD,0,2);

        r.data.F1.(F1Name{ii}).R_std = nanstd(r.data.F1.(F1Name{ii}).All.R,0,2);
        r.data.F1.(F1Name{ii}).Rsum_std = nanstd(r.data.F2.(F2Name{ii}).All.Rsum,0,2);

    end







    %% Plots
    figure(FigNum);clf
    subplot(1,2,1)
    if i1 == 1
        hold on;
    end
    scatter(unique(r.data.PlotParam),r.data.N(:,:),'filled');
    plot_format(ParamName,'Number','',12);
    sgtitle({['{\bf\fontsize{14}' Title '}'],SubTitle});
    grid on
    ylim([0,Inf]);
    legend('F2','F1')

    subplot(1,2,2)
    if i1 == 1
        hold on;
    end
    %         scatter(unique(r.data.PlotParam),r.data.RSum(:,:),100,'ColorVariable',['r','b']);
    %         hold on
    %         ax = scatter(unique(r.data.PlotParam),r.data.R(:,:),20,'filled','ColorVariable',['r','b']);    plot_format(ParamName,'Population','',12);
    ax = scatter(unique(r.data.PlotParam),r.data.RSum(:,:),100,'ColorVariable',['b','r'],'LineWidth',2);
    hold on   
    errorbar(unique(r.data.PlotParam),r.data.R(:,1),r.data.R_std(:,1),'LineStyle','none','Color','b')
    errorbar(unique(r.data.PlotParam),r.data.R(:,2),r.data.R_std(:,2),'LineStyle','none','Color','r')
    grid on;
    plot_format(ParamName,'Population','',12);


    maxlength = max(numel(img(1).clouds),numel(img(2).clouds));
    figure(FigNum+1);clf
    sgtitle({['{\bf\fontsize{14}' Title '}'],SubTitle});
    for ii = 1:numel(img(1).clouds)
        subplot(maxlength,2,ii*2-1);
        if i1 == 1
            hold on;
        end
        scatter(unique(r.data.PlotParam), r.data.F2.(F2Name{ii}).Rsum,'filled')
        hold on
        scatter(unique(r.data.PlotParam),r.data.F2.(F2Name{ii}).R,'filled')
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
        scatter(unique(r.data.PlotParam),r.data.F1.(F1Name{ii}).Rsum,'filled')
        hold on
        scatter(unique(r.data.PlotParam),r.data.F1.(F1Name{ii}).R,'filled')
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