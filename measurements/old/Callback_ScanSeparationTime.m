function Callback_ScanSeparationTime(r)

if r.isInit()
    r.data.param = 15e-3:5e-3:50e-3;
    r.c.setup('var',r.data.param);
    r.data.nlf = nonlinfit.empty;
elseif r.isSet()
    
    r.make('detuning',25.5,'tof',730e-3,'dipole',1.45,'power',0.15,...
        'T',5e-3,'camera','drop 4','asym',250e-6,'tsep',r.data.param(r.c(1)));
    r.upload;
    fprintf(1,'Run %d/%d, Param = %.3f\n',r.c.now,r.c.total,r.data.param(r.c(1))*1e3);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.5);
    
    [nlf,N] = FMI_Analysis;
    r.data.nlf(i1,1) = nlf;
    Ndim = ceil(sqrt(r.c.total));
    figure(10);
    subplot(Ndim,Ndim,i1);
    plot(nlf.x,nlf.y,'.-');
    xlim([3,5]);ylim([0,0.8]);
    title(sprintf('Tsep = %.3f',r.data.param(i1)*1e3));
    
    

end