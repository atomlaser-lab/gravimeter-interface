function Callback_OptimizeFMIFrequency(r)

if r.isInit()
    %Initialize run
    r.data.param = const.randomize(20:0.25:25);
    r.c.setup('var',r.data.param);
    
    r.makerCallback = @makeSequence;
elseif r.isSet()
    %Build/upload/run sequence
    r.make(r.devices.opt.set('detuning',r.data.param(r.c(1))));
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,...
        r.data.param(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    
    [nlf,N] = FMI_Analysis;
    r.data.amp(i1,1) = nlf.c(1,1);
    r.data.contrast(i1,1) = nlf.c(5,1);
    
    figure(10);clf;
    subplot(1,2,1);
    errorbar(r.data.param(1:i1),r.data.amp(1:i1,1),0.025*r.data.amp(1:i1,1),'o');
    plot_format('Detuning [MHz]','Amplitude','',12);
    grid on;
    subplot(1,2,2);
    plot(r.data.param(1:i1),r.data.contrast,'o');
    plot_format('Detuning [MHz]','Contrast','',12);
    grid on;

end