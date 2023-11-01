function Callback_MeasureRamanFreq(r)

if r.isInit()
    r.data.freq = const.randomize([-500:5:500]);
    
    r.c.setup('var',r.data.freq);
elseif r.isSet()
    k = 2*pi*384229441689483/const.c;
    r.make(r.devices.opt.set('raman_df',2*k*9.795*50e-3/(2*pi*1e6) + r.data.freq(r.c(1))*1e-3));
    r.upload;
    fprintf(1,'Run %d/%d, F = %.3f kHz\n',r.c.now,r.c.total,...
        r.data.freq(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);
    img = Abs_Analysis('last');
    if ~img.raw.status.ok()
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
    plot(r.data.freq(1:i1),r.data.N(1:i1,:),'o');
    plot_format('Freq [kHz]','Number','',12);
%     h = legend('m = -1','m = 0','m = 1');
%     set(h,'Location','West');
    title(' Raman frequency using fit over OD')
    grid on
    hold on;
    ylim([0,2e6]);
    
    subplot(1,2,2)
    plot(r.data.freq(1:i1),r.data.R(1:i1,:),'o');
    hold off;
    plot_format('Freq [kHz]','Population','',12);
% %     h = legend('m = -1','m = 0','m = 1');
% %     set(h,'Location','West');
%     title(' Raman frequency using ROI')
    grid on;
%     if r.c.done
%         tNow = datestr(now);
%         caption = sprintf('Determination of Raman frequency %s', tNow);
%         sgtitle(caption)
%     end
end