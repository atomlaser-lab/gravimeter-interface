function Callback_MeasureRamanPulseDuration(r)

if r.isInit()
    r.data.duration = 0:10:150;
    
    r.c.setup('var',r.data.duration);
elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.duration(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, T = %.0f us\n',r.c.now,r.c.total,...
        r.data.duration(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5 + 0.5*rand);
    img = Abs_Analysis('last');
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif i1 > 1 && strcmpi(img.raw.files.name,r.data.files{i1 - 1}.name)
        r.c.decrement;
        return;
    end
    %
    % Store raw data
    %
    r.data.files{i1,1} = img.raw.files;
    %
    % Get processed data
    %
    r.data.N(i1,:) = img.get('N');
    r.data.Nsum(i1,:) = img.get('Nsum');
    r.data.R(i1,:) = r.data.N(i1,:)./sum(r.data.N(i1,:));
    r.data.Rsum(i1,:) = r.data.Nsum(i1,:)./sum(r.data.Nsum(i1,:));
    
    figure(97);clf;
    subplot(1,2,1)
    plot(r.data.duration(1:i1),r.data.N(1:i1,:),'o');
    plot_format('Pulse duration [us]','Number','',12);
%     h = legend('m = -1','m = 0','m = 1');
%     set(h,'Location','West');
    title(' Raman frequency using fit over OD')
    grid on
    hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    plot(r.data.duration(1:i1),r.data.Rsum(1:i1,:),'o');
    hold off;
    plot_format('Pulse duration [us]','Population','',12);
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
        nlf = nonlinfit(r.data.duration(1:r.c(1)),r.data.Rsum(:,1),1e-2);
        nlf.setFitFunc(@(A,R,D,x) A*(1 - 4*R.^2./(4*R^2+D.^2).*sin(2*pi*sqrt(4*R^2+D.^2).*x/2).^2));
        nlf.bounds2('A',[0.7,1,0.95],'R',[0,10,0.1]*1e-3,'D',[-0.1,0.1,0]);
        nlf.fit;
        fprintf(1,'Rabi frequency: %.3f kHz, Detuning = %.3f kHz\n',nlf.c(2,1)*1e3,nlf.c(3,1)*1e3);
        subplot(1,2,2);
        hold on
        xplot = linspace(min(nlf.x),max(nlf.x),1e2);
        plot(xplot,nlf.f(xplot),'-');
        hold off;
        
    end
end