function Callback_Rhys_2D_av_1State(r)

% % Inputs
ClearImage = 1;
FigNum = 5;
Title = 'Pumping';
SubTitle = '';

Param = 1:2:5;
PlotParam = Param;
ParamName = ScanableParameters.TwoPhoton;
Param2 = 0:1:2;
Averages = 4;

if r.isInit()
    %
    % Define parameters
    %
    r.data.param1 = const.randomize(Param);
    r.data.param2 = Param2;
    % Number of averages to do
    r.data.num_avgs = Averages;
    %
    % When to average controls when averaging is done.  If it is 'first',
    % then num_avgs repetitions are done of each set of parameters (param1,
    % param2, etc.) before moving onto the next set.  If it is 'last', then
    % all parameters are scanned first before being repeated num_avgs
    % times.
    %
    r.data.when_to_average = 'first';

    if isavgfirst(r)
        r.c.setup([r.data.num_avgs,numel(r.data.param1),numel(r.data.param2)]);
    else
        r.c.setup([numel(r.data.param1),numel(r.data.param2),r.data.num_avgs]);
    end

    if ClearImage == 1
        figure(FigNum);clf
        figure(FigNum+1);clf
    end

elseif r.isSet()
    if isavgfirst(r)
        i1 = r.c(2);
        i2 = r.c(3);
    else
        i1 = r.c(1);
        i2 = r.c(2);
    end
    r.make(r.devices.opt,'params',r.data.param1(i1),r.data.param2(i2));
    r.upload;
    fprintf(1,'Run %d/%d (%s), Param1 = %.3f, Param2 = %.3f\n',...
        r.c.now,r.c.total,r.c.print,r.data.param1(i1),r.data.param2(i2));
    
elseif r.isAnalyze()
    if isavgfirst(r)
        i1 = r.c(2);
        i2 = r.c(3);
        iavg = r.c(1);
    else
        i1 = r.c(1);
        i2 = r.c(2);
        iavg = r.c(3);
    end
    pause(0.5 + 0.5*rand);

    img = Abs_Analysis_DualState_RT('last');

    if r.c.now == 1
        % Create a structure for each ROI for each image
        F2field_names = arrayfun(@(x) sprintf('ROI%d', x), 1:numel(img(1).clouds), 'UniformOutput', false);
        F1field_names = arrayfun(@(x) sprintf('ROI%d', x), 1:numel(img(2).clouds), 'UniformOutput', false);
        r.data.F2 = cell2struct(cell(size(F2field_names)), F2field_names, 2);
        r.data.F1 = cell2struct(cell(size(F1field_names)), F1field_names, 2);
    end

    if ~img(1).raw.status.ok()
        % 
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence        
        %
        r.c.decrement;
        return;
    elseif r.c.now > 1 && strcmpi(img(1).raw.files.name,r.data.files{r.c.now - 1}.name)
        r.c.decrement;
        return;
    end
    %
    % Store raw data
    %
    r.data.files{i1,i2,iavg} = img(1).raw.files;

    %
    % Get processed data
    %
    r.data.raw.N(i1,i2,iavg) = img(1).clouds(F2_ROI).N;
    r.data.raw.Nsum(i1,i2,iavg) = img(1).clouds(F2_ROI).Nsum;


    %
    % Compute averaging as necessary
    %
    if isavgfirst(r) && r.c.done(1)
        names = fieldnames(r.data.raw)';
        for nn = 1:numel(names)
            p = names{nn};
            r.data.(p).mean(i1,i2,:) = squeeze(mean(r.data.raw.(p)(i1,i2,:,:),3));
            r.data.(p).err(i1,i2,:) = squeeze(std(r.data.raw.(p)(i1,i2,:,:),0,3))./sqrt(r.data.num_avgs);
        end
    elseif ~isavgfirst(r) && r.c.done(2)
        names = fieldnames(r.data.raw)';
        for nn = 1:numel(names)
            p = names{nn};
            r.data.(p).mean = squeeze(mean(r.data.raw.(p),3));
            r.data.(p).err = squeeze(std(r.data.raw.(p),0,3))./sqrt(iavg);
        end
    elseif ~isavgfirst(r)
        names = fieldnames(r.data.raw)';
        for nn = 1:numel(names)
            p = names{nn};
            if iavg > 1
                n1_end = r.c.final(1);
                n2_end = r.c.final(2);
            else
                n1_end = i1;
                n2_end = i2;
            end
            for n1 = 1:n1_end
                for n2 = 1:n2_end
                    if n1 <= i1 && n2 <= i2
                        r.data.(p).mean(n1,n2,:) = squeeze(mean(r.data.raw.(p)(n1,n2,1:iavg,:),3));
                        r.data.(p).err(n1,n2,:) = squeeze(std(r.data.raw.(p)(n1,n2,1:iavg,:),0,3))/sqrt(iavg);
                    else
                        r.data.(p).mean(n1,n2,:) = squeeze(mean(r.data.raw.(p)(n1,n2,1:(iavg - 1),:),3));
                        r.data.(p).err(n1,n2,:) = squeeze(std(r.data.raw.(p)(n1,n2,1:(iavg - 1),:),0,3))/sqrt(iavg - 1);
                    end
                end
            end
        end
    end

    %
    % Plot data
    %
    figure(131);
    if isavgfirst(r) && r.c.done(1)
        clf;
        subplot(1,2,1);
        for nn = 1:size(r.data.N.mean,3)
            h = errorbar(r.data.param1(1:i1),r.data.N.mean(1:i1,i2,nn),r.data.N.err(1:i1,i2,nn),'o');
            h.MarkerFaceColor = h.Color;
            hold on
        end
        xlabel('Param1');ylabel('N');
        
        subplot(1,2,2);
        for nn = 1:size(r.data.R.mean,3)
            h = errorbar(r.data.param1(1:i1),r.data.R.mean(1:i1,i2,nn),r.data.R.err(1:i1,i2,nn),'o');
            h.MarkerFaceColor = h.Color;
            hold on
        end
        xlabel('Param1');ylabel('Population');
    elseif ~isavgfirst(r)
        clf;
        subplot(1,2,1);
        for nn = 1:size(r.data.N.mean,3)
            h = errorbar(r.data.param1(1:i1),r.data.N.mean(1:i1,i2,nn),r.data.N.err(1:i1,i2,nn),'o');
            h.MarkerFaceColor = h.Color;
            hold on
        end
        xlabel('Param1');ylabel('N');
        
        subplot(1,2,2);
        for nn = 1:size(r.data.R.mean,3)
            h = errorbar(r.data.param1(1:i1),r.data.R.mean(1:i1,i2,nn),r.data.R.err(1:i1,i2,nn),'o');
            h.MarkerFaceColor = h.Color;
            hold on
        end
        xlabel('Param1');ylabel('Population');
    end
    
end

end

function res = isavgfirst(r)
    res = strcmpi(r.data.when_to_average,'first');
end