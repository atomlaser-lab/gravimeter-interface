function Callback_MeasureGradiometerFringes(r)

if r.isInit()
%     r.data.run = 1:150;
    r.c.setup(Inf);
elseif r.isSet()
    
%     r.make(0,216.6e-3,1.5,0.215,r.data.phase(r.c(1)),r.data.T(r.c(2)));
    r.make(r.devices.opt.set('phase',mod(r.c(1),360)));
    r.upload;
    fprintf(1,'Run %d/%d\n',r.c.now,r.c.total);
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);

    [out1,out2,~,~,dout] = Gradiometer_Analysis_Symmetric_FMI;
    r.data.N1(i1,:) = out1.N;
    r.data.N2(i1,:) = out2.N;
    r.data.R1(i1,1) = out1.R;
    r.data.R2(i1,1) = out2.R;

    r.data.d(i1,1) = dout;
    
    figure(97);clf;
    subplot(1,2,1);
    plot(1:i1,r.data.R1,'o-',1:i1,r.data.R2,'sq-');
    plot_format('Run','Ratios','',12);
    subplot(1,2,2);
    plot(r.data.R1,r.data.R2,'o');
    plot_format('R1','R2','',12);

    
end