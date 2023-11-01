function Callback_MeasureRamanPulseDuration(r)

if r.isInit()
    r.data.pulse = 2e-6:2e-6:60e-6; %want to see 
    
    r.c.setup('var',r.data.pulse);
elseif r.isSet()
    
    k = 2*pi*384229441689483/const.c;
    r.make(r.devices.opt.set('raman_df',2*k*9.795*50e-3/(2*pi*1e6) - 76.6e-3,'raman_width',r.data.pulse(r.c(1))));

    while 1
        try
            r.upload;
            break;
        catch
%             current_try = current_try + 1;
        end
    end
    fprintf(1,'Run %d/%d, T= %.3f us\n',r.c.now,r.c.total,...
        r.data.pulse(r.c(1))*1e6);
    
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
%     subplot(1,2,1)
    plot(r.data.pulse(1:i1)*1e6,r.data.R(1:i1,:),'o');
    plot_format('Duration [us]','Population','',12);
%     h = legend('m = -1','m = 0','m = 1');
%     set(h,'Location','West');
    title(' Raman frequency using fit over OD')
    grid on
    hold on;
%     ylim([0,2e6]);
  
end