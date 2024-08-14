function Callback_OptimizeFMIDrive(r)

if r.isInit()
    %Initialize run
    r.data.param = const.randomize(0.05:0.05:1);
    r.c.setup('var',r.data.param);
    
    r.makerCallback = @makeSequence;
elseif r.isSet()
    %Build/upload/run sequence
    r.make(r.devices.opt);
    r.devices.d.lockin.driveAmp.set(r.data.param(r.c(1))).write;
    r.upload;
    %Print information about current run
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,...
        r.data.param(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    
    [nlf,N] = FMI_Analysis;
    r.data.N(i1,1) = N.N2;
    r.data.amp(i1,1) = nlf.c(4,1);
    
    figure(10);clf;
    subplot(1,2,1);
    errorbar(r.data.param(1:i1),r.data.N(1:i1,1),0.05*r.data.N(1:i1,1),'o');
    plot_format('Drive [V]','Fitted signal','',12);
    grid on;
    subplot(1,2,2);
    errorbar(r.data.param(1:i1),r.data.amp(1:i1,1),0.025*r.data.amp(1:i1,1),'o');
    plot_format('Drive [V]','Amplitude','',12);
    grid on;

end