function GravImagingFreq(r)

if r.isInit()
    r.data.tof = 30e-3;
    r.data.freq = const.randomize(8:0.1:9);

    r.c.setup('var',r.data.freq);
elseif r.isSet()
    
    r.make(r.data.freq(r.c(1)),r.data.tof,2);
    r.upload;
    fprintf(1,'Run %d/%d, Freq: %.3f V\n',r.c.now,r.c.total,r.data.freq(r.c(1)));
    
elseif r.isAnalyze()
    i1 = r.c(1);
    pause(0.1);
    img = Abs_Analysis('last');
    if ~img(1).raw.status.ok()
        %
        % Checks for an error in loading the files (caused by a missed
        % image) and reruns the last sequence
        %
        r.c.decrement;
        return;
    end
    r.data.N(i1,1) = img.get('N');
    r.data.Nsum(i1,1) = img.get('Nsum');
    
    nlf = nonlinfit(r.data.freq(1:i1),r.data.Nsum(1:i1,1)/1e6,0.05*r.data.Nsum(1:i1,1)/1e6+0.2,r.data.Nsum(1:i1,1)>100e6);
    nlf.setFitFunc(@(A,x0,G,y0,x) A./(1+4*((x-x0)/G).^2)+y0);
    nlf.bounds([0,7.5,0,-0.01],[100,9.5,1,0.1],[1,8.5,0.5,0]);
    figure(10);clf;
    if i1 < 5
        errorbar(nlf.x,nlf.y,nlf.dy,'o');
        plot_format('Voltage [V]','Number of atoms \times 10^6','',12);
    else
        nlf.fit;
        nlf.plot;
        r.data.nlf = nlf;
        fprintf(1,'Resonance: %.3f +/- %.3f, Width: %.3f +/- %.3f\n',nlf.c(2,:),nlf.c(3,:));
    end
    
    r.data.nlf = nlf;
    
end