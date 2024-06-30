function Callback_Rhys_Average(r)

% % Inputs
ClearImage = 1;
FigNum = 5;

Title = 'MOT Scan';
SubTitle = 'Mot coil = 5 A';
NumAverages = 20;
Param = [0.5:0.5:5, 5:1:10];


PlotFactor = 1;
ParamName = ScanableParameters.TwoPhoton;
% ParamName = 'EW Bias (V)';
% % % If there are multiple ROIs, what do you want to count?

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
    end

elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.Param(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, T = %.0f us\n',r.c.now,r.c.total,...
        r.data.Param(r.c(1)));

elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis_GUI('last');

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
    % Get processed data for input region of interest
    %
    AvIndex = r.data.AvIndex(r.c.i);
    ParIndex = r.data.ParamIndex(r.c.i);


    if r.c.i == 1
        % Pre-allocate data
        r.data.All.N = nan(numel(Param),NumAverages);
        r.data.All.NSum = nan(numel(Param),NumAverages);
        r.data.All.OD= nan(numel(Param),NumAverages);
        r.data.All.Width = nan(numel(Param),2,NumAverages);
        r.data.All.T = nan(numel(Param),2,NumAverages);
    end

    % All data
    r.data.All.N(ParIndex,AvIndex) = img(1).clouds.N;
    r.data.All.NSum(ParIndex,AvIndex) = img(1).clouds.Nsum;
    r.data.All.OD(ParIndex,AvIndex) = img(1).clouds.peakOD;
    r.data.All.Width(ParIndex,:,AvIndex) = img(1).clouds.gaussWidth;
    r.data.All.T(ParIndex,:,AvIndex) = img(1).clouds.T;
    % Averaged data
    if NumAverages ~= 1
        r.data.N = nanmean(r.data.All.N,2);
        r.data.NSum = nanmean(r.data.All.NSum,2);
        r.data.OD = nanmean(r.data.All.OD,2);
        r.data.Width = nanmean(r.data.All.Width,3);
        r.data.T = nanmean(r.data.All.T,3);
        % Std
        r.data.N_std = nanstd(r.data.All.N,0,2);
        r.data.NSum_std = nanstd(r.data.All.NSum,0,2);
        r.data.OD = nanstd(r.data.All.OD,0,2);
        r.data.Width = nanstd(r.data.All.Width,0,3);
        r.data.T = nanstd(r.data.All.T,0,3);
    else
        r.data.N = r.data.All.N;
        r.data.NSum = r.data.All.NSum;
        r.data.OD = r.data.All.OD;
        r.data.Width = r.data.All.Width;
        r.data.T = r.data.All.T,;
    end
    %
    % Get processed data for all regions of interest

    %% Plots
    figure(FigNum);clf
%     if NumAverages ~= 1
%         errorbar(unique(r.data.PlotParam),r.data.NSum(:,1),r.data.NSum_std(:,1),'-s','MarkerSize',10,'LineStyle','none','Color','r')
%         hold on
%         errorbar(unique(r.data.PlotParam),r.data.N(:,1),r.data.N_std(:,1),'-s','MarkerSize',10,'LineStyle','none','Color','b')
%         
%     else
%         scatter(unique(r.data.PlotParam),r.data.N(:,1),'Color','b')
%         hold on
%         scatter(unique(r.data.PlotParam),r.data.NSum(:,1),'Color','r')
%         
%     end
%     sgtitle({['{\bf\fontsize{14}' Title '}'],SubTitle});
%     grid on

end
end