function Callback_MeasureBraggPower(r)

if r.isInit()
    r.data.power = 0:0.05:0.9;
%     r.data.t0 = [30e-3,50e-3,75e-3,100e-3];
    r.data.t0 = 50e-3;
    
    r.c.setup('var',r.data.power,r.data.t0);
elseif r.isSet()
    
%     r.make(25.5,730e-3,1.48,r.data.power(r.c(1)),0,130e-3);

    r.make(r.devices.opt,'bragg',{'power',r.data.power(r.c(1)),'t0',r.data.t0(r.c(2))});
    r.upload;
    fprintf(1,'Run %d/%d, P = %.3f\n',r.c.now,r.c.total,...
        r.data.power(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.5);
    img = Abs_Analysis('last');
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    elseif r.c.now > 1 && strcmpi(r.data.files{r.c.now - 1}.name,img.raw.files.name)
        pause(10);
        r.c.decrement;
        return;
    end
    %
    % Store raw data
    %
%     r.data.c{i1,1} = c;
    r.data.files{i1,i2} = img.raw.files;
    %
    % Get processed data
    %
    r.data.N(i1,i2,:) = img.get('N');
    r.data.Nsum(i1,i2,:) = img.get('Nsum');
    r.data.R(i1,i2) = r.data.N(i1,i2,1)./sum(r.data.N(i1,i2,:),3);
    r.data.Rsum(i1,i2) = r.data.Nsum(i1,i2,1)./sum(r.data.Nsum(i1,i2,:),3);

%     [~,N] = FMI_Analysis;
%     r.data.N(i1,:) = [N.N1,N.N2];
%     r.data.R(i1,1) = N.R;
    
    figure(98);clf;
    subplot(1,2,1);
    h = plot(r.data.power(1:i1),r.data.R(1:i1,i2),'o');
    for nn = 1:numel(h)
        set(h(nn),'MarkerFaceColor',h(nn).Color);
    end
%     hold off;
    plot_format('Power [arb units]','Population','Bragg power scan at pulse width of 30 us FWHM',12);
    grid on;
%     h = legend('Slow','Fast');
%     h = legend('-2k','0k','2k','4k');
%     set(h,'Location','West');
    ylim([0,1]);
    xlim([0,1]);
    if r.c.done(1)
        subplot(1,2,2);
        cla;
        for nn = 1:i2
            plot(r.data.power(1:i1),r.data.R(1:i1,nn),'o');
            hold on
        end
        plot_format('Power [arb units]','Population','Bragg power scan at pulse width of 30 us FWHM',12);
        grid on;
    end
    
end