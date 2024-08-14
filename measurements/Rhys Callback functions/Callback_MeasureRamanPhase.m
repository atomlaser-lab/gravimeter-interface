function Callback_MeasureRamanPhase(r)
FigNum = 1500;
% XLabel = 'Run Number';
XLabel = 'Phase [deg]';


if r.isInit()
    r.data.phase = const.randomize(repmat(0:20:180,8,1));
    
    r.c.setup('var',r.data.phase);
elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.phase(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, ph = %.0f deg\n',r.c.now,r.c.total,...
        r.data.phase(r.c(1)));
    
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
%     r.data.R2(i1,:) = r.data.N(i1,[3,7])./sum(r.data.N(i1,[3,7]));
    r.data.R2 = r.data.R;
    
    figure(FigNum);
    subplot(1,2,1)
    scatter(r.data.phase(1:i1),r.data.N(1:i1,:),'filled'); %,'o'
    plot_format(XLabel,'Number','',12);
%     h = legend('m = -1','m = 0','m = 1');
%     set(h,'Location','West');
    title(' Raman frequency using fit over OD')
    grid on
%     hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    scatter(r.data.phase(1:i1),r.data.R2(1:i1,:),'filled'); %,'o'
    hold off;
    plot_format(XLabel,'Population','',12);
    ylim([0,1])
% %     h = legend('m = -1','m = 0','m = 1');
% %     set(h,'Location','West');
%     title(' Raman frequency using ROI')
    grid on;
%     if r.c.done
%         tNow = datestr(now);
%         caption = sprintf('Determination of Raman frequency %s', tNow);
%         sgtitle(caption)
%     end
    if r.c(1) > 4
        [ph,P,dP] = average_repeated_runs(r.data.phase(1:r.c(1)),r.data.R2(:,1));
        dP(dP == 0) = 0.02;
        nlf = nonlinfit(ph,P,dP);
%         nlf.ex = sum(r.data.N,2) > 1e7;
        nlf.setFitFunc(@(y0,C,x0,x) y0 + C/2*sind(2*(x - x0)));
        nlf.bounds2('y0',[0,1,0.5],'C',[0,1,0.8],'x0',[-180,180,0]);
        nlf.fit;
        subplot(1,2,2);
        cla;
        errorbar(nlf.x,nlf.y,nlf.dy,'o');
        fill_markers(gca);
        grid on;
        plot_format(XLabel,'Population','',12);
        ylim([0,1])
        hold on
        xplot = linspace(min(nlf.x),max(nlf.x),1e2);
        plot(xplot,nlf.f(xplot),'-');
        hold off;
        r.data.nlf = nlf;
    end
end