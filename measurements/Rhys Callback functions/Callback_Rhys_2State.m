function Callback_Rhys_2State(r)

% % Inputs
ClearImage = 1;
FigNum = 5;
%,\tau = 20 us
Title = 'Pumping: P_C = 1.5*0.4 mW, P_S/P_C = 2, 3x Mag, t_0 = 0 ms, \Delta = 4 GHz, \delta = -20 - 7.8e-3';
SubTitle = '';
% SubTitle = 'Two-Photon Scan';
% SubTitle = 'Pulse Duration Scan';

Param = 0:40:800;
% Param = 0:1:10;
PlotParam = Param;
ParamName = ScanableParameters.PulseDuration;
% ParamName = 'Phase';
% % % If there are multiple ROIs, what do you want to count?
F2_ROI = 3;
F1_ROI = 2;



if r.isInit()
    r.data.Param = Param;
    r.data.PlotParam = PlotParam;
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
    
    %
    % Get processed data for all regions of interest
    %
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







    %% Plots
    figure(FigNum);clf
    subplot(1,2,1)
    if i1 == 1
        hold on;
    end
    scatter(r.data.PlotParam(1:i1),r.data.N(1:i1,:),'filled');
    plot_format(ParamName,'Number','',12);
    sgtitle({['{\bf\fontsize{14}' Title '}'],SubTitle});  
    grid on
    ylim([0,Inf]);
    legend('F2','F1')
    
    subplot(1,2,2)
    if i1 == 1
        hold on;
    end
    scatter(r.data.PlotParam(1:i1),r.data.Rsum(1:i1,:),100,'ColorVariable',['r','b']);
    hold on
    ax = scatter(r.data.PlotParam(1:i1),r.data.R(1:i1,:),20,'filled','ColorVariable',['r','b']);    plot_format(ParamName,'Population','',12);
    grid on;


    maxlength = max(numel(img(1).clouds),numel(img(2).clouds));
    figure(FigNum+1);clf
    sgtitle({['{\bf\fontsize{14}' Title '}'],SubTitle});  
    for ii = 1:numel(img(1).clouds)
        subplot(maxlength,2,ii*2-1);
        if i1 == 1
            hold on;
        end
        scatter(r.data.PlotParam(1:i1), r.data.F2.(F2Name{ii}).Rsum,'filled')
        hold on
        scatter(r.data.PlotParam(1:i1),r.data.F2.(F2Name{ii}).R,'filled')
%         if max(r.data.F2.(F2Name{ii}).Rsum) < 0.1
%             ylim([0,0.1])
%         else
%             ylim([0,1])
%         end
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
%         if max(r.data.F1.(F1Name{ii}).Rsum) < 0.1
%             ylim([0,0.1])
%         else
%             ylim([0,1])
%         end        
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