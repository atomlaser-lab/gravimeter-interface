function GravEvapOptimize(r)

if r.isInit()
    %Initialize run
    r.data.evapEnd = const.randomize(7.65:0.05:7.85);
    r.data.evapRate = 0.1:0.05:0.3;
    r.c.setup('var',r.data.evapEnd,r.data.evapRate);
    
    r.makerCallback = @makeSequence;
elseif r.isSet()
    %Build/upload/run sequence
    r.make(r.devices.opt.set('params',[r.data.evapEnd(r.c(1)),r.data.evapRate(r.c(2))]));
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, Evap End = %.3f V, Evap Rate = %.2f V/s\n',r.c.now,r.c.total,...
        r.data.evapEnd(r.c(1)),r.data.evapRate(r.c(2)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.5); %Wait for other image analysis program to finish with files
    %Analyze image data from last image
    img = Abs_Analysis('last');
    if ~img.raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end
    r.data.files{i1,i2} = img.raw.files;
    r.data.N(i1,i2) = img.get('N');
    r.data.Nsum(i1,i2) = img.get('Nsum');
%     r.data.T(i1,i2) = sqrt(prod(c.T));
%     r.data.OD(i1,i2) = c.peakOD;
    
    figure(10);
    subplot(1,2,1);
    errorbar(r.data.evapEnd(1:i1),r.data.N(1:i1,i2)/1e6,0.025*r.data.N(1:i1,i2)/1e6,'o');
    plot_format('Start voltage','Number of atoms \times 10^6','',12);
    ylim([0,3]);
    grid on;
    if r.c.done(1)
        subplot(1,2,2);
        errorbar(r.data.evapEnd,r.data.N(:,i2)/1e6,0.025*r.data.N(:,i2)/1e6,'o');
        plot_format('Start voltage','Number of atoms \times 10^6','',12);
        hold on;
        r.data.str{i2} = sprintf('Rate = %.2f',r.data.evapRate(i2));
        legend(r.data.str);
        ylim([0,3]);
        grid on;
    end
end