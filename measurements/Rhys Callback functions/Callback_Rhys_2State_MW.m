function Callback_Rhys_2State_MW(r)

% % % Inputs
ClearImage = 0;
FigNum = 5;
TOF = 25e-3;

Title = 'In Trap MW Transfer';
Param = -0.3:0.1:0.3;
ParamName = 'df (kHz)';

% % % If there are multiple ROIs, what do you want to count?
F2_ROI = 3;
F1_ROI = 2;

if r.isInit()
    r.data.freq1 = const.f_Rb_groundHFS - 315e3 - 8e3;
    r.data.freq2 = const.f_Rb_groundHFS + Param*1e3;
    r.data.Param = Param;
    r.c.setup('var',r.data.Param);

elseif r.isSet()
    if numel(r.data.freq1) > 1
        r.devices.mku.writeList(r.data.freq1(r.c(1))/2,r.data.freq2/2);
    elseif numel(r.data.freq2) > 1
        r.devices.mku.writeList(r.data.freq1/2,r.data.freq2(r.c(1))/2);   
    end
    r.make(r.devices.opt, 'tof', TOF).upload;
    fprintf(1,'Run %d/%d, df = %.0f kHz\n',r.c.now,r.c.total,...
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

        r.data.F2.(F2Name{ii}).R(i1,:) = img(1).clouds(ii).N/(sum(vertcat(img(1).clouds(:).N)) + sum(vertcat(img(2).clouds(:).N)));
        r.data.F2.(F2Name{ii}).Rsum(i1,:) = img(1).clouds(ii).Nsum/(sum(vertcat(img(1).clouds(:).Nsum)) + sum(vertcat(img(2).clouds(:).Nsum)));
    end

    F1Name = fieldnames(r.data.F1);
    for ii = 1:size(img(2).clouds,1)
        r.data.F1.(F1Name{ii}).N(i1,:) = img(1).clouds(ii).N;
        r.data.F1.(F1Name{ii}).Nsum(i1,:) = img(1).clouds(ii).Nsum;

        r.data.F1.(F2Name{ii}).R(i1,:) = img(2).clouds(ii).N/(sum(vertcat(img(1).clouds(:).N)) + sum(vertcat(img(2).clouds(:).N)));
        r.data.F1.(F2Name{ii}).Rsum(i1,:) = img(2).clouds(ii).Nsum/(sum(vertcat(img(1).clouds(:).Nsum)) + sum(vertcat(img(2).clouds(:).Nsum)));        
    end

    figure(FigNum);clf
    subplot(1,2,1)
    scatter(r.data.Param(1:i1),r.data.N(1:i1,:),'filled');
    plot_format(ParamName,'Number','',12);
    title(' Raman frequency using fit over OD')
    grid on
    ylim([0,Inf]);
    legend('F2','F1')
    
    subplot(1,2,2)
    ax = scatter(r.data.Param(1:i1),r.data.R(1:i1,:));
    plot_format(ParamName,'Population','',12);
    grid on;


    maxlength = max(numel(img(1).clouds),numel(img(2).clouds));
    figure(FigNum+1);clf
    sgtitle('All ROI')
    for ii = 1:numel(img(1).clouds)
        subplot(maxlength,2,ii*2-1)
        scatter(r.data.Param(1:i1), r.data.F2.(F2Name{ii}).Rsum)
        hold on
        scatter(r.data.Param(1:i1),r.data.F2.(F2Name{ii}).R)
        if max(r.data.F2.(F2Name{ii}).Rsum) < 0.5
            ylim([0,0.5])
        else
            ylim([0,1])
        end
        title(sprintf('ROI %g',ii))
        ylabel('Pop')
        if ii == numel(img(1).clouds)
            xlabel(ParamName)
        end
    end
    for ii = 1:numel(img(2).clouds)
        subplot(maxlength,2,ii*2)
        scatter(r.data.Param(1:i1),r.data.F1.(F1Name{ii}).Rsum)
        hold on
        scatter(r.data.Param(1:i1),r.data.F1.(F1Name{ii}).R)
        if max(r.data.F1.(F1Name{ii}).Rsum) < 0.5
            ylim([0,0.5])
        else
            ylim([0,1])
        end        
        title(sprintf('ROI %g',ii))
        ylabel('Pop')
        if ii == numel(img(1).clouds)
            xlabel(ParamName)
        end
    end

end