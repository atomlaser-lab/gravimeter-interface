function Callback_MeasureRamanFreq(r)

if r.isInit()
    r.data.freq = 0e-3 + 1e-3*const.randomize(-200:1:200);
    
    r.c.setup('var',r.data.freq);
elseif r.isSet()
    r.make(r.devices.opt,'params',r.data.freq(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, F = %.3f MHz\n',r.c.now,r.c.total,...
        r.data.freq(r.c(1)));
    
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
    
    figure(97);clf;
    subplot(1,2,1)
    scatter(r.data.freq(1:i1),r.data.N(1:i1,:),'filled');
    plot_format('Freq [MHz]','Number','',12);
%     h = legend('m = -1','m = 0','m = 1');
%     set(h,'Location','West');
    title(' Raman frequency using fit over OD')
    grid on
    hold on;
    ylim([0,Inf]);
    
    subplot(1,2,2)
    ax = scatter(r.data.freq(1:i1),r.data.R(1:i1,:));
    ax(6).Marker="x";
    ax(7).Marker="x";
    ax(8).Marker="x";
    hold off;
    plot_format('Freq [MHz]','Population','',12);
    h = legend('2,2','2,1,','2,0','2,-1','2,-2','1,-1','1,0','1,1');
    set(h,'Location','West');
%     title(' Raman frequency using ROI')
    grid on;
%     if r.c.done
%         tNow = datestr(now);
%         caption = sprintf('Determination of Raman frequency %s', tNow);
%         sgtitle(caption)
%     end
%     if r.c(1) > 4
%         nlf = nonlinfit(r.data.freq(1:r.c(1)),r.data.R(:,2),1e-2);
%         nlf.setFitFunc(@(A,R,x0,x) A*(1 - 4*R.^2./(4*R^2+(x-x0).^2).*sin(2*pi*sqrt(4*R^2+(x-x0).^2).*100/2).^2));
%         [~,idx] = min(nlf.y);
%         nlf.bounds2('A',[0.5,2,0.95],'R',[0,10,0.3]*1e-3,'x0',[min(nlf.x),max(nlf.x),nlf.x(idx)]);
%         nlf.fit;
%         fprintf(1,'Rabi frequency: %.3f kHz, Center = %.3f kHz\n',nlf.c(2,1)*1e3,nlf.c(3,1)*1e3);
%         subplot(1,2,2);
%         hold on
%         xplot = linspace(min(nlf.x),max(nlf.x),1e2);
%         plot(xplot,nlf.f(xplot),'-');
%         hold off;
%         
%     end
end