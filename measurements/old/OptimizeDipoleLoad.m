function OptimizeDipoleLoad(r)

if r.isInit()
    %Initialize run
    r.data.P25 = const.randomize(5:2.5:15);
    r.data.P50 = 10:5:25;
    r.c.setup('var',r.data.P25,r.data.P50);
    
    r.makerCallback = @makeSequence;
elseif r.isSet()
    %Build/upload/run sequence
    r.make(8.5,35e-3,2,r.data.P25(r.c(1)),r.data.P50(r.c(2)));
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, P25 = %.1f W, P50 = %.1f W\n',r.c.now,r.c.total,...
        r.data.P25(r.c(1)),r.data.P50(r.c(2)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    i2 = r.c(2);
    pause(0.5); %Wait for other image analysis program to finish with files
    %Analyze image data from last image
    c = Abs_Analysis('last');
    r.data.files{i1,i2} = {c.raw.files(1).name,c.raw.files(2).name};
    r.data.N(i1,i2) = c.N;
    r.data.Nsum(i1,i2) = c.Nsum;
    r.data.T(i1,i2) = sqrt(prod(c.T));
    r.data.OD(i1,i2) = c.peakOD;
    r.data.PSD(i1,i2) = c.PSD;
    
    figure(11);
    subplot(1,2,1);
    errorbar(r.data.P25(1:i1),r.data.N(1:i1,i2)/1e6,0.025*r.data.N(1:i1,i2)/1e6,'o');
    plot_format('P25 [W]','Number of atoms \times 10^6','',12);
    ylim([0,5]);
    grid on;
    if i1 == r.c.imax(1)
        subplot(1,2,2);
        errorbar(r.data.P25,r.data.N(:,i2)/1e6,0.025*r.data.N(:,i2)/1e6,'o');
        plot_format('P25 [W]','Number of atoms \times 10^6','',12);
        hold on;
        r.data.str{i2} = sprintf('P50 = %.1f',r.data.P50(i2));
        legend(r.data.str);
        ylim([0,5]);
        grid on;
    end
end