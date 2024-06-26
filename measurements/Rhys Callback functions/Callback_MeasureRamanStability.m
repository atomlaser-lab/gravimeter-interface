function Callback_MeasureRamanStability(r)
FigNum = 99;
XLabel = 'Run Number';
% XLabel = 'Pulse Duration [us]';


if r.isInit()
    r.data.run = 1:100;
    
    r.c.setup('var',r.data.run);
elseif r.isSet()
    r.make(r.devices.opt);
    r.upload;
    fprintf(1,'Run %d/%d\n',r.c.now,r.c.total);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis_DualState_RT('last');
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img(1).raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end
    %
    % Store raw data
    %
    r.data.files{i1,1} = img(1).raw.files;
    %
    % Get processed data
    %
    r.data.N(i1,:) = [img(1).get('N'),img(2).get('N')];
    r.data.Nsum(i1,:) = [img(1).get('Nsum'),img(2).get('Nsum')];
%     if numel(img) > 1
%         r.data.N(i1,2) = r.data.N(i1,2) - r.data.N(i1,1);
%         r.data.Nsum(i1,2) = r.data.Nsum(i1,2) - r.data.Nsum(i1,1);
%     end
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));
    r.data.R2(i1,:) = r.data.N(i1,[3,7])./sum(r.data.N(i1,[3,7]));
    
    figure(FigNum);
    subplot(1,2,1)
    scatter(r.data.run(1:i1),r.data.N(1:i1,:),'filled'); %,'o'
    plot_format(XLabel,'Number','',12);
%     h = legend('m = -1','m = 0','m = 1');
%     set(h,'Location','West');
    title(' Raman frequency using fit over OD')
    grid on
%     hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    scatter(r.data.run(1:i1),r.data.R2(1:i1,:),'filled'); %,'o'
    hold off;
    plot_format(XLabel,'Population','',12);
    ylim([0,1])
% %     h = legend('m = -1','m = 0','m = 1');
% %     set(h,'Location','West');
%     title(' Raman frequency using ROI')
    grid on;

end